import re
import streamlit as st
import pandas as pd
import plotly.graph_objects as go

st.set_page_config(page_title="PCR 프라이머 분석", page_icon="🔬", layout="wide")

# 염기서열 입력에서 공백과 줄바꿈을 제거합니다.
def clean_sequence(value):
    return re.sub(r"\s+", "", value).upper()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]


def gc_percent(sequence):
    return 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)


def wallace_tm(sequence):
    # 간이 Wallace 식: 짧은 프라이머의 대략적인 Tm 추정에만 사용합니다.
    at = sequence.count("A") + sequence.count("T")
    gc = sequence.count("G") + sequence.count("C")
    return 2 * at + 4 * gc


def find_all(sequence, pattern):
    # 겹치는 일치 구간도 모두 찾습니다(0부터 시작하는 인덱스).
    return [m.start() for m in re.finditer(f"(?={re.escape(pattern)})", sequence)]


def validate_sequence(sequence, label):
    if not sequence:
        return f"{label}을(를) 입력해 주세요."
    invalid = sorted(set(sequence) - set("ATGC"))
    if invalid:
        return f"{label}에 허용되지 않는 문자가 있습니다: {', '.join(invalid)} (A/T/G/C만 지원)"
    return None


st.title("🔬 PCR 프라이머 분석")
st.write("DNA 서열과 프라이머를 입력해 **GC 함량, 간이 Tm, 결합 위치 및 예상 PCR 증폭 산물**을 확인합니다.")
st.caption("교육용 분석 도구입니다. 프라이머 결합은 완전 일치만 탐색하며, 실제 PCR 성공을 보장하지 않습니다.")

# 바로 작동을 확인할 수 있도록 가상 예시 서열을 제공합니다.
example_forward = "ATGCGTACGTTAGCTAGTCA"
example_reverse = "CGATCGTACCTAGGCTATGC"
example_template = example_forward + "GACTACGATCGATGCTAGCTGATCGTAGCTAGCTACGATCGATGCTAGCTGATC" + reverse_complement(example_reverse)

with st.expander("예시 데이터 사용 방법", expanded=False):
    st.write("아래 입력창에는 가상의 테스트 데이터가 기본으로 들어 있습니다. 그대로 분석 버튼을 누르면 예상 증폭 산물을 확인할 수 있습니다.")

with st.form("pcr_analysis_form"):
    template_input = st.text_area("표적 DNA 서열 (5′ → 3′)", value=example_template, height=140,
                                  help="공백과 줄바꿈은 자동 제거됩니다. A/T/G/C만 지원합니다.")
    left, right = st.columns(2)
    with left:
        forward_input = st.text_input("정방향 프라이머 (5′ → 3′)", value=example_forward)
    with right:
        reverse_input = st.text_input("역방향 프라이머 (5′ → 3′)", value=example_reverse,
                                      help="역방향 프라이머도 실제 주문하는 서열 방향인 5′ → 3′으로 입력합니다.")
    analyze = st.form_submit_button("PCR 분석하기", type="primary", use_container_width=True)

if analyze:
    template = clean_sequence(template_input)
    forward = clean_sequence(forward_input)
    reverse = clean_sequence(reverse_input)

    errors = [validate_sequence(s, label) for s, label in [
        (template, "표적 DNA"), (forward, "정방향 프라이머"), (reverse, "역방향 프라이머")
    ]]
    errors = [e for e in errors if e]
    if len(template) > 20000:
        errors.append("표적 DNA는 한 번에 20,000 bp 이하로 입력해 주세요.")
    if any(len(s) > 60 for s in (forward, reverse)):
        errors.append("프라이머 길이는 60 nt 이하로 입력해 주세요.")

    if errors:
        for error in errors:
            st.error(error)
    else:
        st.subheader("1. 프라이머 기본 특성")
        primer_info = pd.DataFrame([
            {"프라이머": "정방향", "길이 (nt)": len(forward), "GC 함량 (%)": round(gc_percent(forward), 2),
             "간이 Tm (°C)": wallace_tm(forward)},
            {"프라이머": "역방향", "길이 (nt)": len(reverse), "GC 함량 (%)": round(gc_percent(reverse), 2),
             "간이 Tm (°C)": wallace_tm(reverse)},
        ])
        st.dataframe(primer_info, hide_index=True, use_container_width=True)
        st.caption("Tm은 Wallace 식(2×A/T 수 + 4×G/C 수)의 단순 추정치입니다. 특히 긴 프라이머나 실제 반응 조건에서 정확한 Tm으로 사용하면 안 됩니다.")

        warnings = []
        for name, primer in [("정방향", forward), ("역방향", reverse)]:
            if not 18 <= len(primer) <= 25:
                warnings.append(f"{name} 프라이머 길이가 일반적인 초기 검토 범위(18~25 nt)를 벗어납니다.")
            if not 40 <= gc_percent(primer) <= 60:
                warnings.append(f"{name} 프라이머 GC 함량이 일반적인 초기 검토 범위(40~60%)를 벗어납니다.")
            if not 14 <= len(primer) <= 25:
                warnings.append(f"{name} 프라이머의 Wallace Tm은 해당 길이에서 특히 부정확할 수 있습니다.")
        if abs(wallace_tm(forward) - wallace_tm(reverse)) > 5:
            warnings.append("두 프라이머의 간이 Tm 차이가 5°C를 초과합니다.")
        for warning in warnings:
            st.warning(warning)

        st.subheader("2. 결합 위치 및 증폭 산물 예측")
        # 정방향 프라이머는 입력된 가닥과 같은 서열을 찾습니다.
        # 역방향 프라이머는 입력된 5′→3′ 서열의 역상보 서열을 찾습니다.
        forward_sites = find_all(template, forward)
        reverse_binding_sequence = reverse_complement(reverse)
        reverse_sites = find_all(template, reverse_binding_sequence)

        col1, col2, col3 = st.columns(3)
        col1.metric("DNA 길이", f"{len(template):,} bp")
        col2.metric("정방향 결합 위치 수", len(forward_sites))
        col3.metric("역방향 결합 위치 수", len(reverse_sites))

        # 정방향 결합 위치가 역방향 결합 위치보다 앞에 있어야
        # 두 프라이머가 서로 마주 보는 방향으로 증폭 산물을 만듭니다.
        amplicons = []
        for f_start in forward_sites:
            for r_start in reverse_sites:
                if f_start + len(forward) <= r_start:
                    amplicons.append({
                        "시작 위치 (1-based)": f_start + 1,
                        "끝 위치 (1-based)": r_start + len(reverse),
                        "예상 길이 (bp)": r_start + len(reverse) - f_start,
                        "_f_start": f_start,
                        "_r_start": r_start,
                    })

        if not forward_sites or not reverse_sites:
            st.warning("한쪽 또는 양쪽 프라이머의 완전 일치 결합 위치를 찾지 못했습니다.")
        elif not amplicons:
            st.warning("각 프라이머의 일치 위치는 있지만, 서로 마주 보는 방향의 증폭 산물은 찾지 못했습니다.")
        else:
            st.success(f"예상 PCR 증폭 산물 {len(amplicons)}개를 찾았습니다.")
            st.dataframe(pd.DataFrame(amplicons).drop(columns=["_f_start", "_r_start"]),
                         hide_index=True, use_container_width=True)
            if len(amplicons) > 1:
                st.warning("여러 증폭 산물이 예측됩니다. 표적 서열 안에서 결합 위치가 중복되는지 확인하세요.")

        # 염기서열 위치를 선분으로 표시하여 프라이머와 산물의 위치를 비교합니다.
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[1, len(template)], y=[0, 0], mode="lines",
                                 line=dict(color="gray", width=8), name="표적 DNA"))
        for i, pos in enumerate(forward_sites):
            fig.add_trace(go.Scatter(x=[pos + 1, pos + len(forward)], y=[0.18, 0.18],
                                     mode="lines", line=dict(color="#16a34a", width=12),
                                     name="정방향 프라이머" if i == 0 else None,
                                     showlegend=(i == 0)))
        for i, pos in enumerate(reverse_sites):
            fig.add_trace(go.Scatter(x=[pos + 1, pos + len(reverse)], y=[-0.18, -0.18],
                                     mode="lines", line=dict(color="#ea580c", width=12),
                                     name="역방향 프라이머 결합 부위" if i == 0 else None,
                                     showlegend=(i == 0)))
        for i, amp in enumerate(amplicons[:20]):
            fig.add_trace(go.Scatter(x=[amp["시작 위치 (1-based)"], amp["끝 위치 (1-based)"]],
                                     y=[0.4 + 0.09 * i, 0.4 + 0.09 * i], mode="lines",
                                     line=dict(color="#2563eb", width=5),
                                     name="예상 증폭 산물" if i == 0 else None,
                                     showlegend=(i == 0),
                                     hovertemplate=f"예상 길이: {amp['예상 길이 (bp)']} bp<extra></extra>"))
        fig.update_layout(title="프라이머 결합 위치와 예상 증폭 구간", xaxis_title="DNA 위치 (bp)",
                          yaxis=dict(visible=False), height=350, margin=dict(l=20, r=20, t=60, b=30))
        st.plotly_chart(fig, use_container_width=True)
        if len(amplicons) > 20:
            st.caption("그래프에는 앞의 20개 증폭 산물만 표시했습니다. 표에는 전체 결과가 나옵니다.")

        st.subheader("3. 결과 해석 및 한계")
        if len(amplicons) == 1:
            st.info(f"입력한 표적 서열에서 완전 일치 기준 예상 증폭 산물은 {amplicons[0]['예상 길이 (bp)']} bp입니다.")
        elif len(amplicons) == 0:
            st.info("입력한 표적 서열에서 완전 일치 기준 예상 증폭 산물을 찾지 못했습니다.")
        st.write("이 분석은 입력된 **한 개의 DNA 서열**에서 완전 일치 결합만 찾습니다. "
                 "프라이머 이합체·헤어핀, 불일치 결합, 전체 유전체에서의 비특이적 증폭, "
                 "염 농도·반응 온도 등은 평가하지 않습니다. 따라서 실제 실험의 증폭 성공 여부를 판단할 수 없습니다.")
