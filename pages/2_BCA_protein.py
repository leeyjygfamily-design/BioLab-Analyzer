import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import io
import json
import zipfile
import xml.etree.ElementTree as ET

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
            # 질량 환산은 선택 기능입니다. 세 항목을 모두 입력했을 때만 계산합니다.
            fields = ["전체 시료 질량 (g)", "추출에 사용한 시료 질량 (g)", "회수 상등액 부피 (mL)"]
            entries = [row.get(field, np.nan) for field in fields]
            present = [not (pd.isna(v) or str(v).strip() == "") for v in entries]
            if any(present) and not all(present):
                raise ValueError(f"{label} {idx + 1}행: 전체 질량 환산을 하려면 질량 2개와 상등액 부피를 모두 입력해주세요.")
            if all(present):
                nums = [pd.to_numeric(v, errors="coerce") for v in entries]
                if any(pd.isna(v) or not np.isfinite(v) or v <= 0 for v in nums):
                    raise ValueError(f"{label} {idx + 1}행: 질량과 상등액 부피는 0보다 큰 숫자여야 합니다.")
                if nums[1] > nums[0]:
                    raise ValueError(f"{label} {idx + 1}행: 추출에 사용한 질량이 전체 시료 질량보다 클 수 없습니다.")
                result.update(dict(zip(fields, map(float, nums))))
            else:
                result.update({field: np.nan for field in fields})
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
        # 원액 농도는 BCA 측정 전에 희석한 배수를 반영한 추출액 농도입니다.
        # 추출된 양 = 추출액 농도 × 회수한 상등액의 전체 부피
        # 전체 시료 환산 = 추출된 양 × (전체 시료 질량 / 추출에 사용한 질량)
        if np.isfinite(row["회수 상등액 부피 (mL)"]):
            protein_mg = row["원액 농도 (mg/mL)"] * row["회수 상등액 부피 (mL)"]
            total_mg = protein_mg * (row["전체 시료 질량 (g)"] / row["추출에 사용한 시료 질량 (g)"])
            row["추출된 단백질 양 (mg)"] = protein_mg if measured >= 0 else np.nan
            row["전체 시료 추정 단백질 양 (mg)"] = total_mg if measured >= 0 else np.nan
            row["전체 시료 추정 단백질 비율 (%)"] = (total_mg / (row["전체 시료 질량 (g)"] * 1000) * 100) if measured >= 0 else np.nan
        else:
            for field in ("추출된 단백질 양 (mg)", "전체 시료 추정 단백질 양 (mg)", "전체 시료 추정 단백질 비율 (%)"):
                row[field] = np.nan
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


# -------------------- SMBL 데이터 불러오기 --------------------
# SMBL 파일의 저장 형식은 버전/측정 모드에 따라 달라질 수 있습니다.
# 명확한 파장-흡광도 쌍이 있는 구조화 데이터만 자동 인식합니다.
WAVELENGTH_MIN, WAVELENGTH_MAX = 561.0, 563.0


def _number(value):
    try:
        v = float(str(value).strip())
        return v if np.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def _field_type(name):
    name = str(name).lower().replace(" ", "").replace("_", "")
    if any(t in name for t in ("wavelength", "lambda", "파장", "nm")):
        return "wavelength"
    if any(t in name for t in ("absorbance", "absorb", "흡광도", "opticaldensity")) or name in ("abs", "a"):
        return "absorbance"
    return None


def _extract_rows(rows, source):
    """표 헤더가 파장/흡광도임을 명확히 나타낼 때만 추출."""
    out = []
    if not rows:
        return out
    for row in rows:
        if not isinstance(row, dict):
            continue
        wave_keys = [k for k in row if _field_type(k) == "wavelength"]
        abs_keys = [k for k in row if _field_type(k) == "absorbance"]
        for wk in wave_keys:
            w = _number(row[wk])
            if w is None:
                continue
            for ak in abs_keys:
                a = _number(row[ak])
                if a is not None:
                    out.append({"파일/데이터": source, "파장 (nm)": w, "흡광도": a})
    return out


def _extract_json(obj, source, depth=0):
    if depth > 15:
        return []
    result = []
    if isinstance(obj, list):
        result.extend(_extract_rows(obj, source))
        for value in obj:
            if isinstance(value, (dict, list)):
                result.extend(_extract_json(value, source, depth + 1))
    elif isinstance(obj, dict):
        result.extend(_extract_rows([obj], source))
        wave_arrays = [(k, v) for k, v in obj.items() if _field_type(k) == "wavelength" and isinstance(v, list)]
        abs_arrays = [(k, v) for k, v in obj.items() if _field_type(k) == "absorbance" and isinstance(v, list)]
        for _, ws in wave_arrays:
            for _, av in abs_arrays:
                if len(ws) == len(av):
                    result.extend(_extract_rows([{"wavelength": w, "absorbance": a} for w, a in zip(ws, av)], source))
        for value in obj.values():
            if isinstance(value, (dict, list)):
                result.extend(_extract_json(value, source, depth + 1))
    return result


def _extract_xml(data, source):
    """Vernier SMBL의 DataSet별 파장/흡광도 열을 대응시켜 읽습니다."""
    root = ET.fromstring(data)
    result = []
    for dataset in root.findall('.//DataSet'):
        dataset_name = (dataset.findtext('DataSetName') or '이름 없는 데이터').strip()
        wave_values, absorb_values = None, None
        for column in dataset.findall('.//DataColumn'):
            kind = _field_type(column.findtext('DataObjectName') or '')
            raw = (column.findtext('ColumnCells') or '').strip()
            if not raw:
                continue
            try:
                values = [float(v) for v in raw.split()]
            except ValueError:
                continue
            if kind == 'wavelength':
                wave_values = values
            elif kind == 'absorbance':
                absorb_values = values
        if wave_values is None or absorb_values is None or len(wave_values) != len(absorb_values):
            continue
        # 각 DataSet은 하나의 측정 스펙트럼이므로 가장 가까운 파장 하나만 선택합니다.
        candidates = [(abs(w - 562.0), w, a) for w, a in zip(wave_values, absorb_values)
                      if np.isfinite(w) and np.isfinite(a) and WAVELENGTH_MIN <= w <= WAVELENGTH_MAX]
        if candidates:
            _, wavelength, absorbance = min(candidates, key=lambda v: v[0])
            result.append({'파일/데이터': source, '데이터 세트': dataset_name,
                           '파장 (nm)': wavelength, '흡광도': absorbance})
    if result:
        return result
    # 다른 XML 파일에서는 기존의 행 기반 인식 방식을 사용합니다.
    rows = []
    for element in root.iter():
        row = {k.split('}')[-1]: v for k, v in element.attrib.items()}
        for child in element:
            if len(child) == 0 and child.text:
                row[child.tag.split('}')[-1]] = child.text
        if row:
            rows.append(row)
    return _extract_rows(rows, source)


def _extract_table(data, source):
    try:
        text = data.decode('utf-8-sig')
    except UnicodeDecodeError:
        return []
    for sep in (None, '\t', ',', ';'):
        try:
            df = pd.read_csv(io.StringIO(text), sep=sep, engine='python', on_bad_lines='skip')
            if len(df) and len(df.columns) >= 2:
                found = _extract_rows(df.to_dict('records'), source)
                if found:
                    return found
        except (ValueError, pd.errors.ParserError, UnicodeError):
            continue
    return []


def _extract_file(data, name):
    """ZIP, JSON, XML, 표 형식으로 저장된 측정 파일을 탐색합니다."""
    if len(data) > 25 * 1024 * 1024:
        raise ValueError('파일은 25 MB 이하만 지원합니다.')
    parts = [(name, data)]
    if zipfile.is_zipfile(io.BytesIO(data)):
        parts = []
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for item in z.infolist():
                if item.is_dir() or item.file_size > 10 * 1024 * 1024:
                    continue
                if item.filename.lower().endswith(('.json', '.xml', '.csv', '.tsv', '.txt')):
                    parts.append((f'{name} / {item.filename}', z.read(item)))
    result = []
    for part_name, part in parts:
        stripped = part.lstrip()[:1]
        try:
            if stripped in (b'{', b'['):
                result.extend(_extract_json(json.loads(part.decode('utf-8-sig')), part_name))
            elif stripped == b'<':
                result.extend(_extract_xml(part, part_name))
            else:
                result.extend(_extract_table(part, part_name))
        except (ValueError, ET.ParseError, UnicodeDecodeError, RecursionError):
            continue
    # XML/JSON 중첩 탐색 등으로 중복 추출한 완전히 동일한 레코드를 제거
    return pd.DataFrame(result).drop_duplicates().reset_index(drop=True) if result else pd.DataFrame(
        columns=['파일/데이터', '파장 (nm)', '흡광도'])


def _apply_import(df, selected, kind, row_number, measurement_col):
    target = df.copy().reset_index(drop=True)
    while len(target) <= row_number:
        blank = {c: np.nan for c in target.columns}
        if kind == '미지 시료':
            blank.update({'시료명': f'미지 시료 {len(target) + 1}', '희석배수': 1.0})
        else:
            blank['단백질 농도 (mg/mL)'] = np.nan
        target = pd.concat([target, pd.DataFrame([blank])], ignore_index=True)
    target.loc[row_number, measurement_col] = float(selected)
    return target


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
- **전체 시료 환산:** 원액 농도 × 회수 상등액 전체 부피 × (전체 시료 질량 ÷ 추출에 사용한 질량)입니다. BCA 반응액의 부피가 아닌 상등액 회수량을 사용합니다.
- **주의:** CV 10%는 이 앱에서 사용하는 *점검용 예시 기준*이며, 실제 허용 기준은 실험법·장비·시료에 따라 다릅니다. R²가 높아도 결과의 정확성이 보장되지는 않습니다.
- **SMBL 업로드:** 파일 구조가 자동 인식되는 경우 561~563 nm 범위의 흡광도만 선택할 수 있습니다. 측정 전 파장을 562 nm로 설정하세요.
- **예시 데이터:** 아래 초기값은 사용 방법을 보여주기 위한 가상 데이터입니다. 실제 실험값으로 교체해주세요.
""")

# 표 데이터는 session_state에 보관하여 SMBL 값 적용 후에도 편집 상태를 유지합니다.
if 'bca_standard_table' not in st.session_state:
    st.session_state.bca_standard_table = pd.DataFrame({
        '단백질 농도 (mg/mL)': [0.0, 0.25, 0.5, 1.0, 1.5, 2.0],
        '흡광도 1': [0.00, 0.10, 0.21, 0.41, 0.62, 0.82],
        '흡광도 2': [np.nan] * 6, '흡광도 3': [np.nan] * 6,
    })
if 'bca_unknown_table' not in st.session_state:
    st.session_state.bca_unknown_table = pd.DataFrame({
        '시료명': ['미지 시료 A'], '흡광도 1': [0.30],
        '흡광도 2': [np.nan], '흡광도 3': [np.nan], '희석배수': [1.0],
        '전체 시료 질량 (g)': [np.nan], '추출에 사용한 시료 질량 (g)': [np.nan],
        '회수 상등액 부피 (mL)': [np.nan],
    })
if 'bca_editor_version' not in st.session_state:
    st.session_state.bca_editor_version = 0

st.subheader('① 분광광도계 파일 불러오기 (선택)')
st.info('BCA 측정 시 분광광도계의 파장을 **562 nm**로 설정하세요. 파일 분석에서는 **561~563 nm (562 ± 1 nm)**를 허용합니다.')
st.caption('Vernier SpectroVis Plus 측정 파일(.smbl)을 업로드해 보세요. CSV/TSV/JSON/XML 형식도 지원합니다. 파일 형식에 따라 자동 추출이 불가능할 수 있습니다.')
uploaded = st.file_uploader('측정 파일 업로드', type=['smbl', 'csv', 'tsv', 'json', 'xml', 'txt'], accept_multiple_files=True)
if uploaded:
    signature = tuple((f.name, f.size, f.getvalue()[:128]) for f in uploaded)
    if st.session_state.get('bca_upload_signature') != signature:
        frames, failed = [], []
        for f in uploaded:
            try:
                extracted = _extract_file(f.getvalue(), f.name)
                if extracted.empty:
                    failed.append(f.name)
                else:
                    frames.append(extracted)
            except (ValueError, OSError, zipfile.BadZipFile) as exc:
                failed.append(f'{f.name} ({exc})')
        st.session_state.bca_import_data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        st.session_state.bca_upload_signature = signature
        st.session_state.bca_import_failed = failed
    for f in st.session_state.get('bca_import_failed', []):
        st.warning(f'{f}: 파장-흡광도 데이터 구조를 자동 인식하지 못했습니다. CSV로 내보내거나 실제 SMBL 파일 구조를 확인해야 합니다.')
    extracted = st.session_state.bca_import_data
    if not extracted.empty:
        valid = extracted[extracted['파장 (nm)'].between(WAVELENGTH_MIN, WAVELENGTH_MAX)].copy().reset_index(drop=True)
        st.caption(f'파장-흡광도 쌍 {len(extracted)}개 인식 / 허용 범위(561~563 nm) {len(valid)}개')
        if valid.empty:
            st.error('561~563 nm 범위의 측정값이 없습니다. 측정 파장과 파일 내용을 확인하세요.')
        else:
            st.success(f'562 ± 1 nm 범위의 흡광도 {len(valid)}개를 찾았습니다.')
            st.dataframe(valid, hide_index=True, use_container_width=True)
            st.caption('각 데이터 세트에서 562 nm와 가장 가까운 파장 하나를 선택했습니다. 데이터 세트 이름과 실제 시료를 확인하세요.')
            # 표준용액과 미지 시료를 각각 여러 개 일괄 입력할 수 있습니다.
            if '데이터 세트' in valid.columns:
                with st.expander('SMBL 데이터 세트 일괄 입력', expanded=True):
                    st.write('파일 안의 모든 데이터 세트를 표준용액 또는 미지 시료 입력표로 가져옵니다.')
                    batch_kind = st.radio('일괄 입력 대상', ['표준용액', '미지 시료'],
                                          horizontal=True, key='bca_batch_kind')
                    if batch_kind == '표준용액':
                        st.caption('데이터 세트 이름이 실제 농도인 경우에만 선택하세요. 문자 이름은 농도로 변환하지 않습니다.')
                        unit = st.selectbox('데이터 세트 이름의 농도 단위',
                                            ['µg/mL (mg/mL로 변환)', 'mg/mL (변환 없음)'])
                    else:
                        st.caption('데이터 세트 이름을 그대로 시료명으로 사용합니다. 숫자 이름도 시료명으로 사용할 수 있습니다.')
                        dilution = st.number_input('일괄 적용할 희석배수 (이후 개별 수정 가능)',
                                                    min_value=1.0, value=1.0, step=1.0)
                    if st.button('모든 데이터 세트를 입력표에 적용', type='secondary'):
                        rows = []
                        invalid = []
                        grouped = {}
                        # 같은 데이터 세트 이름이 반복되면 최대 3회 반복 측정값으로 묶습니다.
                        for _, item in valid.iterrows():
                            name = str(item['데이터 세트']).strip()
                            if not name or name.lower() == 'nan':
                                name = str(item['파일/데이터']).split('/')[-1].strip()
                            a = float(item['흡광도'])
                            if not np.isfinite(a) or a < 0:
                                invalid.append(f'{name}: 음수 또는 잘못된 흡광도')
                                continue
                            if batch_kind == '표준용액':
                                conc = _number(name)
                                if conc is None or conc < 0:
                                    invalid.append(f'{name}: 농도 이름이 숫자가 아님')
                                    continue
                                if unit.startswith('µg'):
                                    conc /= 1000
                                key = conc
                            else:
                                key = name
                            grouped.setdefault(key, []).append(a)
                        if any(len(values) > 3 for values in grouped.values()):
                            st.error('같은 이름으로 측정한 값이 3개를 초과합니다. 현재 입력표는 시료당 최대 3회 측정만 지원합니다.')
                        elif invalid:
                            st.error('일괄 입력할 수 없는 데이터가 있습니다: ' + '; '.join(invalid))
                        elif not grouped:
                            st.error('일괄 입력할 측정값이 없습니다.')
                        else:
                            for key, values in grouped.items():
                                row = {'흡광도 1': values[0],
                                       '흡광도 2': values[1] if len(values) > 1 else np.nan,
                                       '흡광도 3': values[2] if len(values) > 2 else np.nan}
                                if batch_kind == '표준용액':
                                    row['단백질 농도 (mg/mL)'] = float(key)
                                else:
                                    row['시료명'] = str(key)
                                    row['희석배수'] = float(dilution)
                                    row.update({'전체 시료 질량 (g)': np.nan, '추출에 사용한 시료 질량 (g)': np.nan, '회수 상등액 부피 (mL)': np.nan})
                                rows.append(row)
                            if batch_kind == '표준용액':
                                key = 'bca_standard_table'
                                result = pd.DataFrame(rows).sort_values('단백질 농도 (mg/mL)').reset_index(drop=True)
                            else:
                                key = 'bca_unknown_table'
                                result = pd.DataFrame(rows)
                            st.session_state[key] = result
                            st.session_state.bca_editor_version += 1
                            st.session_state.pop('bca_result', None)
                            st.rerun()
            with st.form('bca_import_form'):
                choices = [f"{i+1}. {row['파일/데이터']} | {row['파장 (nm)']:.2f} nm | 흡광도 {row['흡광도']:.5f}" for i, row in valid.iterrows()]
                chosen = st.selectbox('입력할 측정값', range(len(choices)), format_func=lambda i: choices[i])
                target_type = st.radio('입력 대상', ['표준용액', '미지 시료'], horizontal=True)
                current_df = st.session_state.bca_standard_table if target_type == '표준용액' else st.session_state.bca_unknown_table
                row_index = st.number_input('입력할 행 번호 (1부터 시작)', min_value=1, max_value=100, value=1, step=1)
                column = st.selectbox('흡광도 입력 칸', ['흡광도 1', '흡광도 2', '흡광도 3'])
                submitted = st.form_submit_button('선택한 흡광도 입력표에 적용')
            if submitted:
                value = valid.iloc[chosen]['흡광도']
                if value < 0:
                    st.error('음수 흡광도는 현재 BCA 분석 입력 규칙에서 허용하지 않습니다.')
                else:
                    key = 'bca_standard_table' if target_type == '표준용액' else 'bca_unknown_table'
                    st.session_state[key] = _apply_import(st.session_state[key], value, target_type, int(row_index)-1, column)
                    if target_type == '미지 시료':
                        item = valid.iloc[chosen]
                        dataset_name = str(item.get('데이터 세트', '')).strip()
                        if dataset_name and dataset_name.lower() != 'nan':
                            st.session_state[key].loc[int(row_index)-1, '시료명'] = dataset_name
                    st.session_state.bca_editor_version += 1
                    st.session_state.pop('bca_result', None)
                    st.rerun()
else:
    st.session_state.pop('bca_upload_signature', None)
    st.session_state.pop('bca_import_data', None)
    st.session_state.pop('bca_import_failed', None)

st.subheader('② 표준용액 데이터')
st.caption('같은 농도에서 최대 3회 측정한 흡광도를 한 행에 입력하세요. 측정하지 않은 칸은 비워두면 됩니다.')
standard_data = st.data_editor(st.session_state.bca_standard_table, num_rows='dynamic', use_container_width=True,
                               hide_index=True, key=f"bca_standards_{st.session_state.bca_editor_version}")
st.session_state.bca_standard_table = standard_data.copy()

st.subheader('③ 미지 시료 데이터')
st.caption('희석하지 않은 시료는 희석배수 1을 입력하세요. 전체 단백질 양까지 계산하려면 시료별로 질량 2개와 원심분리 후 회수한 상등액의 전체 부피를 입력하세요. (BCA 반응에 넣은 부피가 아닙니다.)')
# 이전 버전에서 생성된 session_state 입력표에도 새 열을 추가합니다.
for field in ('전체 시료 질량 (g)', '추출에 사용한 시료 질량 (g)', '회수 상등액 부피 (mL)'):
    if field not in st.session_state.bca_unknown_table.columns:
        st.session_state.bca_unknown_table[field] = np.nan
unknown_data = st.data_editor(st.session_state.bca_unknown_table, num_rows='dynamic', use_container_width=True,
                              hide_index=True, key=f"bca_unknowns_{st.session_state.bca_editor_version}")
st.session_state.bca_unknown_table = unknown_data.copy()

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
    st.subheader("④ 분석 결과")
    a, b, c = st.columns(3)
    a.metric("표준곡선 R²", f"{r['r2']:.4f}")
    b.metric("표준용액 수", f"{len(standards)}개")
    c.metric("미지 시료 수", f"{len(unknowns)}개")
    st.code(f"흡광도 = {r['slope']:.5f} × 농도 + {r['intercept']:.5f}", language="text")

    st.subheader("표준곡선과 반복 측정 오차막대")
    fig = go.Figure()
    # 미지 시료별 색상을 사용자가 직접 선택합니다. 그래프 두 곳에 동일하게 적용합니다.
    vivid_colors = {
        '선명한 빨강': '#F5222D',
        '오렌지': '#FF7800',
        '골드': '#E5AD00',
        '초록': '#16A34A',
        '파랑': '#1769FF',
        '보라': '#9333EA',
        '핫핑크': '#FF1493',
        '진한 회색': '#374151',
    }
    color_names = list(vivid_colors)
    st.markdown('**미지 시료 그래프 색상 선택**')
    st.caption('각 시료의 색상을 선택하면 표준곡선과 농도 비교 그래프에 동시에 반영됩니다.')
    selected_colors = []
    color_cols = st.columns(min(3, max(1, len(unknowns))))
    for i, unknown in enumerate(unknowns):
        with color_cols[i % len(color_cols)]:
            selected_name = st.selectbox(
                f"{unknown['시료명']} 색상",
                options=color_names,
                index=i % len(color_names),
                key=f"bca_color_{i}_{unknown['시료명']}",
            )
            selected_colors.append(vivid_colors[selected_name])
    symbols = ['diamond', 'triangle-up', 'square', 'cross', 'star', 'hexagon', 'x', 'pentagon']
    # 표준용액은 모두 같은 하늘색 원으로 표시하고 범례도 하나로 묶습니다.
    # 미지 시료만 서로 다른 색상과 모양으로 구분합니다.
    fig.add_trace(go.Scatter(
        x=[standard['단백질 농도 (mg/mL)'] for standard in standards],
        y=[standard['평균 흡광도'] for standard in standards],
        mode='markers', name='표준용액',
        marker=dict(size=11, color='#87CEEB', symbol='circle'),
        error_y=dict(
            type='data',
            array=[standard['표준편차'] if np.isfinite(standard['표준편차']) else 0 for standard in standards],
            visible=True
        ),
        hovertemplate='농도 %{x:.3f} mg/mL<br>평균 흡광도 %{y:.4f}<extra>표준용액</extra>'
    ))
    grid = np.linspace(r['min'], r['max'], 100)
    fig.add_trace(go.Scatter(x=grid, y=r['slope'] * grid + r['intercept'],
                             mode='lines', name='선형 회귀선', line=dict(color='#444444', dash='dash')))
    for i, unknown in enumerate(unknowns):
        fig.add_trace(go.Scatter(
            x=[unknown['측정 시료 농도 (mg/mL)']], y=[unknown['평균 흡광도']],
            mode='markers', name=unknown['시료명'],
            marker=dict(symbol=symbols[(i+2) % len(symbols)], size=13,
                        color=selected_colors[i], line=dict(color='black', width=1)),
            error_y=dict(type='data', array=[unknown['표준편차'] if np.isfinite(unknown['표준편차']) else 0], visible=True),
            hovertemplate=f"{unknown['시료명']}<br>농도 %{{x:.3f}} mg/mL<br>흡광도 %{{y:.4f}}<extra></extra>"
        ))
    fig.update_layout(xaxis_title="단백질 농도 (mg/mL)", yaxis_title="흡광도", height=420)
    st.plotly_chart(fig, use_container_width=True)
    st.caption("오차막대는 반복 측정의 ±1 표본 표준편차입니다. 반복 측정이 1회이면 오차막대를 표시하지 않습니다.")

    st.subheader("⑤ 표준용액 신뢰도 분석")
    standard_cols = ["단백질 농도 (mg/mL)", "측정 횟수", "평균 흡광도", "표준편차", "CV (%)", "예측 흡광도", "잔차"]
    st.dataframe(table(standards, standard_cols).round(4), use_container_width=True, hide_index=True)
    residual_fig = px.bar(table(standards, standard_cols), x="단백질 농도 (mg/mL)", y="잔차",
                          title="표준곡선 잔차 (실제 평균 흡광도 − 예측 흡광도)")
    residual_fig.add_hline(y=0, line_dash="dash")
    st.plotly_chart(residual_fig, use_container_width=True)

    st.subheader("⑥ 미지 시료별 단백질 농도")
    unknown_cols = ["시료명", "측정 횟수", "평균 흡광도", "표준편차", "CV (%)",
                    "희석배수", "측정 시료 농도 (mg/mL)", "원액 농도 (mg/mL)", "해석"]
    st.dataframe(table(unknowns, unknown_cols).round(4), use_container_width=True, hide_index=True)
    concentration_fig = px.bar(table(unknowns, unknown_cols), x="시료명", y="원액 농도 (mg/mL)",
                               title="미지 시료별 원액 단백질 농도", text_auto=".3f")
    concentration_fig.update_traces(marker_color=selected_colors)
    st.plotly_chart(concentration_fig, use_container_width=True)

    st.subheader("⑦ 전체 시료의 추정 단백질 양")
    st.caption("회수 상등액의 전체 부피를 기준으로 계산합니다. 추출 과정의 단백질 손실은 보정하지 않았으며, 일부 시료의 결과를 전체로 환산하려면 시료가 균일하다는 가정이 필요합니다.")
    mass_cols = ["시료명", "전체 시료 질량 (g)", "추출에 사용한 시료 질량 (g)",
                 "회수 상등액 부피 (mL)", "추출된 단백질 양 (mg)",
                 "전체 시료 추정 단백질 양 (mg)", "전체 시료 추정 단백질 비율 (%)"]
    mass_rows = [u for u in unknowns if np.isfinite(u["전체 시료 추정 단백질 양 (mg)"])]
    if mass_rows:
        st.dataframe(table(mass_rows, mass_cols).round(4), hide_index=True, use_container_width=True)
        mass_fig = px.bar(table(mass_rows, mass_cols), x="시료명", y="전체 시료 추정 단백질 양 (mg)",
                          title="전체 시료별 추정 단백질 양", text_auto=".3f")
        color_lookup = {u["시료명"]: selected_colors[i] for i, u in enumerate(unknowns)}
        mass_fig.update_traces(marker_color=[color_lookup.get(u["시료명"], "#1769FF") for u in mass_rows])
        st.plotly_chart(mass_fig, use_container_width=True)
        if any("표준곡선 농도 범위 밖" in u["해석"] for u in mass_rows):
            st.warning("일부 시료는 표준곡선 범위 밖에서 농도가 계산되었습니다. 전체 단백질 양도 신뢰하기 어려우므로 희석 후 재측정을 고려하세요.")
    else:
        st.info("미지 시료 입력표에 전체 질량, 추출에 사용한 질량, 회수 상등액 부피를 입력하면 계산 결과가 표시됩니다.")

    st.subheader("⑧ 결과 자동 해석")
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
