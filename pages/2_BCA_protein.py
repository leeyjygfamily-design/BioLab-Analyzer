import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(page_title="BCA 단백질 정량", page_icon="🧪", layout="wide")

# -------------------- 분석 함수 --------------------
def measurement_stats(values):
    """반복 측정의 평균, 표본 표준편차, 변동계수를 계산합니다."""
    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1)) if len(arr) >= 2 else np.nan
    cv = float(sd / abs(mean) * 100) if len(arr) >= 2 and abs(mean) > 1e-12 else np.nan
    return mean, sd, cv


def calculate_regression(concentrations, absorbances):
    """기존 기능 유지: 최소제곱 직선의 기울기, 절편, 결정계수."""
    x = np.asarray(concentrations, dtype=float)
    y = np.asarray(absorbances, dtype=float)
    if len(x) < 2 or len(np.unique(x)) < 2:
        raise ValueError("서로 다른 농도의 표준용액이 최소 2개 필요합니다.")
    if np.ptp(y) <= 1e-12:
        raise ValueError("흡광도 변화가 없어 유효한 표준곡선을 만들 수 없습니다.")
    slope, intercept = np.polyfit(x, y, 1)
    if not np.isfinite(slope) or abs(slope) <= 1e-12:
        raise ValueError("회귀선의 기울기가 0에 가까워 농도를 계산할 수 없습니다.")
    predicted = slope * x + intercept
    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(slope), float(intercept), float(r_squared)


def calculate_unknown(absorbance, slope, intercept):
    """기존 기능 유지: 회귀식에서 미지 시료의 측정 농도 계산."""
    if abs(slope) <= 1e-12:
        return None
    return (absorbance - intercept) / slope


def parse_measurements(df, label, has_concentration=False):
    """빈 칸은 무시하되, 숫자가 아닌 값이나 음수는 오류로 처리합니다."""
    rows = []
    for idx, row in df.iterrows():
        name = str(row.get("시료명", "")).strip() if not has_concentration else f"표준용액 {idx + 1}"
        if not has_concentration and (not name or name.lower() == "nan"):
            raise ValueError(f"{label} {idx + 1}행의 시료명을 입력해주세요.")
        concentration = None
        if has_concentration:
            concentration = pd.to_numeric(row.get("단백질 농도 (mg/mL)"), errors="coerce")
            if pd.isna(concentration) or not np.isfinite(concentration) or concentration < 0:
                raise ValueError(f"{label} {idx + 1}행의 농도는 0 이상의 숫자여야 합니다.")
        measurements = []
        for col in ("흡광도 1", "흡광도 2", "흡광도 3"):
            raw = row.get(col)
            if pd.isna(raw) or str(raw).strip() == "":
                continue
            number = pd.to_numeric(raw, errors="coerce")
            if pd.isna(number) or not np.isfinite(number) or number < 0:
                raise ValueError(f"{label} {idx + 1}행의 {col}에는 0 이상의 숫자를 입력해주세요.")
            measurements.append(float(number))
        if not measurements:
            raise ValueError(f"{label} {idx + 1}행에 흡광도를 최소 1개 입력해주세요.")
        mean, sd, cv = measurement_stats(measurements)
        result = {"시료명": name, "측정 횟수": len(measurements), "평균 흡광도": mean,
                  "표준편차": sd, "CV (%)": cv, "측정값": measurements}
        if has_concentration:
            result["단백질 농도 (mg/mL)"] = float(concentration)
        else:
            dilution = pd.to_numeric(row.get("희석배수"), errors="coerce")
            if pd.isna(dilution) or not np.isfinite(dilution) or dilution < 1:
                raise ValueError(f"{label} {idx + 1}행의 희석배수는 1 이상이어야 합니다.")
            result["희석배수"] = float(dilution)
        rows.append(result)
    return rows


def analyze(standards_df, unknowns_df):
    standards = parse_measurements(standards_df, "표준용액", True)
    unknowns = parse_measurements(unknowns_df, "미지 시료", False)
    if len(standards) < 2:
        raise ValueError("표준용액은 최소 2개 필요합니다.")
    concentrations = [r["단백질 농도 (mg/mL)"] for r in standards]
    if len(set(concentrations)) != len(concentrations):
        raise ValueError("표준용액 농도가 중복되었습니다. 같은 농도의 반복 측정값은 한 행에 입력해주세요.")
    if not unknowns:
        raise ValueError("미지 시료를 최소 1개 입력해주세요.")

    # 각 농도의 반복 측정 평균을 이용해 표준곡선을 구합니다.
    # 반복 횟수가 달라도 각 표준 농도에 동일한 가중치를 부여합니다.
    x = np.asarray(concentrations)
    y = np.asarray([r["평균 흡광도"] for r in standards])
    slope, intercept, r2 = calculate_regression(x, y)
    if slope <= 0:
        raise ValueError("표준곡선 기울기가 양수가 아닙니다. 입력 데이터와 실험 조건을 확인해주세요.")
    fitted = slope * x + intercept
    for row, predicted in zip(standards, fitted):
        row["예측 흡광도"] = float(predicted)
        row["잔차"] = row["평균 흡광도"] - float(predicted)

    low, high = float(min(x)), float(max(x))
    for row in unknowns:
        measured = calculate_unknown(row["평균 흡광도"], slope, intercept)
        row["측정 시료 농도 (mg/mL)"] = measured
        row["원액 농도 (mg/mL)"] = measured * row["희석배수"]
        warnings = []
        if measured < low or measured > high:
            warnings.append("표준곡선 농도 범위 밖(외삽)")
        if measured < 0:
            warnings.append("계산 농도가 음수: 정량 불가")
        if np.isfinite(row["CV (%)"]) and row["CV (%)"] > 10:
            warnings.append("반복 측정 CV 10% 초과: 측정값 점검")
        if row["측정 횟수"] == 1:
            warnings.append("반복 측정 없음: 편차 평가 불가")
        row["해석"] = " / ".join(warnings) if warnings else "특이 경고 없음"

    return {"standards": standards, "unknowns": unknowns, "slope": slope,
            "intercept": intercept, "r2": r2, "min": low, "max": high}


def table(rows, columns):
    return pd.DataFrame([{col: row.get(col, np.nan) for col in columns} for row in rows])


# -------------------- 화면 --------------------
st.title("🧪 BCA 단백질 정량 분석")
st.write("표준곡선으로 단백질 농도를 계산하고, 반복 측정의 일관성과 결과 해석의 한계를 확인합니다.")

with st.expander("분석 방법 및 주의사항", expanded=False):
    st.markdown("""
- **표준곡선:** 단백질 농도(x)와 평균 흡광도(y)의 관계를 직선 `y = ax + b`로 나타냅니다.
- **표준편차(SD):** 반복 측정값의 퍼짐 정도입니다. 2회 이상 측정했을 때 표본 표준편차를 계산합니다.
- **CV(%):** 표준편차 ÷ 평균 × 100입니다. 평균이 0이면 계산하지 않습니다.
- **잔차:** 실제 평균 흡광도 − 회귀식 예측 흡광도입니다. 특정 농도에서 직선과 얼마나 차이가 나는지 보여줍니다.
- **원액 농도:** 회귀식으로 구한 측정 시료 농도 × 희석배수입니다.
- **주의:** CV 10%는 이 앱에서 사용하는 *점검용 예시 기준*이며, 실제 허용 기준은 실험법·장비·시료에 따라 다릅니다. R²가 높아도 결과의 정확성이 보장되지는 않습니다.
- **예시 데이터:** 아래 초기값은 사용 방법을 보여주기 위한 가상 데이터입니다. 실제 실험값으로 교체해주세요.
""")

st.subheader("① 표준용액 데이터")
st.caption("같은 농도에서 최대 3회 측정한 흡광도를 한 행에 입력하세요. 측정하지 않은 칸은 비워두면 됩니다.")
standard_default = pd.DataFrame({
    "단백질 농도 (mg/mL)": [0.0, 0.25, 0.5, 1.0, 1.5, 2.0],
    "흡광도 1": [0.00, 0.10, 0.21, 0.41, 0.62, 0.82],
    "흡광도 2": [np.nan] * 6,
    "흡광도 3": [np.nan] * 6,
})
standard_data = st.data_editor(standard_default, num_rows="dynamic", use_container_width=True,
                               hide_index=True, key="bca_standards")

st.subheader("② 미지 시료 데이터")
st.caption("희석하지 않은 시료는 희석배수 1을 입력하세요. 여러 시료를 행으로 추가할 수 있습니다.")
unknown_default = pd.DataFrame({
    "시료명": ["미지 시료 A"], "흡광도 1": [0.30],
    "흡광도 2": [np.nan], "흡광도 3": [np.nan], "희석배수": [1.0]
})
unknown_data = st.data_editor(unknown_default, num_rows="dynamic", use_container_width=True,
                              hide_index=True, key="bca_unknowns")

if st.button("BCA 분석하기", type="primary", use_container_width=True):
    try:
        st.session_state.bca_result = analyze(standard_data, unknown_data)
    except (ValueError, TypeError, np.linalg.LinAlgError) as exc:
        st.session_state.pop("bca_result", None)
        st.error(str(exc))

if "bca_result" in st.session_state:
    r = st.session_state.bca_result
    standards, unknowns = r["standards"], r["unknowns"]
    st.divider()
    st.subheader("③ 분석 결과")
    a, b, c = st.columns(3)
    a.metric("표준곡선 R²", f"{r['r2']:.4f}")
    b.metric("표준용액 수", f"{len(standards)}개")
    c.metric("미지 시료 수", f"{len(unknowns)}개")
    st.code(f"흡광도 = {r['slope']:.5f} × 농도 + {r['intercept']:.5f}", language="text")

    st.subheader("표준곡선과 반복 측정 오차막대")
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[s["단백질 농도 (mg/mL)"] for s in standards],
        y=[s["평균 흡광도"] for s in standards], mode="markers",
        name="표준용액 평균", marker=dict(size=10),
        error_y=dict(type="data", array=[s["표준편차"] if np.isfinite(s["표준편차"]) else 0 for s in standards], visible=True),
        hovertemplate="농도 %{x:.3f} mg/mL<br>평균 흡광도 %{y:.4f}<extra></extra>"
    ))
    grid = np.linspace(r["min"], r["max"], 100)
    fig.add_trace(go.Scatter(x=grid, y=r["slope"] * grid + r["intercept"],
                             mode="lines", name="선형 회귀선"))
    fig.add_trace(go.Scatter(x=[u["측정 시료 농도 (mg/mL)"] for u in unknowns],
                             y=[u["평균 흡광도"] for u in unknowns],
                             mode="markers", name="미지 시료", marker=dict(symbol="diamond", size=10),
                             text=[u["시료명"] for u in unknowns],
                             hovertemplate="%{text}<br>농도 %{x:.3f} mg/mL<br>흡광도 %{y:.4f}<extra></extra>"))
    fig.update_layout(xaxis_title="단백질 농도 (mg/mL)", yaxis_title="흡광도", height=420)
    st.plotly_chart(fig, use_container_width=True)
    st.caption("오차막대는 반복 측정의 ±1 표본 표준편차입니다. 반복 측정이 1회이면 오차막대를 표시하지 않습니다.")

    st.subheader("④ 표준용액 신뢰도 분석")
    standard_cols = ["단백질 농도 (mg/mL)", "측정 횟수", "평균 흡광도", "표준편차", "CV (%)", "예측 흡광도", "잔차"]
    st.dataframe(table(standards, standard_cols).round(4), use_container_width=True, hide_index=True)
    residual_fig = px.bar(table(standards, standard_cols), x="단백질 농도 (mg/mL)", y="잔차",
                          title="표준곡선 잔차 (실제 평균 흡광도 − 예측 흡광도)")
    residual_fig.add_hline(y=0, line_dash="dash")
    st.plotly_chart(residual_fig, use_container_width=True)

    st.subheader("⑤ 미지 시료별 단백질 농도")
    unknown_cols = ["시료명", "측정 횟수", "평균 흡광도", "표준편차", "CV (%)",
                    "희석배수", "측정 시료 농도 (mg/mL)", "원액 농도 (mg/mL)", "해석"]
    st.dataframe(table(unknowns, unknown_cols).round(4), use_container_width=True, hide_index=True)
    concentration_fig = px.bar(table(unknowns, unknown_cols), x="시료명", y="원액 농도 (mg/mL)",
                               title="미지 시료별 원액 단백질 농도", text_auto=".3f")
    st.plotly_chart(concentration_fig, use_container_width=True)

    st.subheader("⑥ 결과 자동 해석")
    st.write(f"표준곡선 농도 범위: **{r['min']:.3f} ~ {r['max']:.3f} mg/mL**")
    if r["r2"] < 0.99:
        st.warning("표준곡선 R²가 0.99 미만입니다. 직선 적합도와 잔차를 점검하세요. 0.99는 참고용 점검 기준입니다.")
    else:
        st.info("표준곡선 R²가 0.99 이상입니다. 다만 이것만으로 정량 결과의 정확성을 보장하지는 않습니다.")
    for u in unknowns:
        if u["해석"] == "특이 경고 없음":
            st.success(f"{u['시료명']}: 표시된 점검 기준에서 특이 경고가 없습니다.")
        else:
            st.warning(f"{u['시료명']}: {u['해석']}")
    if any(s["측정 횟수"] == 1 for s in standards):
        st.info("일부 표준용액은 1회만 측정되어 표준편차·CV를 계산할 수 없습니다.")
    st.caption("이 결과는 교육용 데이터 분석입니다. 실제 실험의 정량 가능 범위, 시약 간섭, 블랭크 처리 및 검량선 검증은 별도로 확인해야 합니다.")
