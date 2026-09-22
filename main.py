import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px


# ---------------------------------------------------------
# 기본 설정
# ---------------------------------------------------------

st.set_page_config(
    page_title="BioLab Analyzer",
    page_icon="🧬",
    layout="wide"
)


# ---------------------------------------------------------
# DNA 분석 함수
# ---------------------------------------------------------

def dna_validation(sequence):
    """DNA 서열에 A, T, G, C 이외의 문자가 있는지 확인합니다."""
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    valid_bases = set("ATGC")
    invalid_bases = set(sequence) - valid_bases

    return sequence, invalid_bases


def count_bases(sequence):
    """DNA 서열의 각 염기 개수를 계산합니다."""
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }


def calculate_gc_content(sequence):
    """DNA 서열의 GC 함량을 계산합니다."""
    if len(sequence) == 0:
        return 0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def complementary_sequence(sequence):
    """DNA의 상보적 염기서열을 만듭니다."""
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    return "".join(complement[base] for base in sequence)


def find_sequence(sequence, target):
    """DNA 서열에서 특정 염기서열이 시작되는 위치를 찾습니다."""
    target = target.upper().replace(" ", "")

    positions = []
    start = 0

    while True:
        position = sequence.find(target, start)

        if position == -1:
            break

        positions.append(position + 1)
        start = position + 1

    return positions


def compare_sequences(sequence1, sequence2):
    """두 DNA 서열의 차이를 비교합니다."""
    max_length = max(len(sequence1), len(sequence2))

    differences = []

    for i in range(max_length):
        base1 = sequence1[i] if i < len(sequence1) else "-"
        base2 = sequence2[i] if i < len(sequence2) else "-"

        if base1 != base2:
            differences.append({
                "위치": i + 1,
                "서열 1": base1,
                "서열 2": base2
            })

    return differences


# ---------------------------------------------------------
# BCA 분석 함수
# ---------------------------------------------------------

def calculate_regression(concentrations, absorbances):
    """표준용액 데이터로 1차 선형회귀식을 계산합니다."""
    x = np.array(concentrations, dtype=float)
    y = np.array(absorbances, dtype=float)

    slope, intercept = np.polyfit(x, y, 1)

    predicted = slope * x + intercept

    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    if ss_tot == 0:
        r_squared = 0
    else:
        r_squared = 1 - (ss_res / ss_tot)

    return slope, intercept, r_squared


def calculate_unknown_concentration(absorbance, slope, intercept):
    """표준곡선의 회귀식을 이용하여 미지 시료 농도를 계산합니다."""
    if slope == 0:
        return None

    concentration = (absorbance - intercept) / slope

    return concentration


# ---------------------------------------------------------
# 세션 상태
# ---------------------------------------------------------

if "page" not in st.session_state:
    st.session_state.page = "home"


# ---------------------------------------------------------
# 공통 함수
# ---------------------------------------------------------

def go_to(page):
    """화면을 이동합니다."""
    st.session_state.page = page


def home_button():
    """홈 화면으로 이동하는 버튼입니다."""
    if st.button("← 홈으로 돌아가기"):
        go_to("home")
        st.rerun()


# ---------------------------------------------------------
# 홈 화면
# ---------------------------------------------------------

def show_home():
    st.title("🧬 BioLab Analyzer")
    st.subheader("생체분자 데이터 분석 프로그램")

    st.write(
        "생명공학 실험에서 얻은 DNA 서열과 "
        "BCA 단백질 정량 데이터를 분석할 수 있는 프로그램입니다."
    )

    st.divider()

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🧬 DNA 서열 분석")
        st.write(
            "DNA 서열의 염기 조성, GC 함량, "
            "상보적 서열 및 서열 차이를 분석합니다."
        )

        if st.button("DNA 분석 시작 →", use_container_width=True):
            go_to("dna_input")
            st.rerun()

    with col2:
        st.subheader("🧪 BCA 단백질 정량")
        st.write(
            "BCA 실험의 표준용액 흡광도 데이터를 이용하여 "
            "미지 시료의 단백질 농도를 계산합니다."
        )

        if st.button("BCA 분석 시작 →", use_container_width=True):
            go_to("bca_input")
            st.rerun()


# ---------------------------------------------------------
# DNA 입력 화면
# ---------------------------------------------------------

def show_dna_input():
    st.title("🧬 DNA 서열 분석")

    home_button()

    st.write("분석할 DNA 서열을 입력하세요.")

    sequence = st.text_area(
        "DNA 서열",
        placeholder="예: ATGCGTACCGTA",
        height=150
    )

    target = st.text_input(
        "찾고 싶은 염기서열",
        placeholder="예: ATG"
    )

    st.write("두 DNA 서열을 비교하려면 아래에 비교할 서열도 입력하세요.")

    sequence2 = st.text_area(
        "비교할 DNA 서열",
        placeholder="예: ATGCGTATCGTA",
        height=120
    )

    if st.button("DNA 분석하기", type="primary", use_container_width=True):

        sequence, invalid_bases = dna_validation(sequence)

        if not sequence:
            st.error("DNA 서열을 입력해주세요.")
            return

        if invalid_bases:
            invalid_text = ", ".join(sorted(invalid_bases))
            st.error(
                f"잘못된 문자가 포함되어 있습니다: {invalid_text}\n\n"
                "DNA 서열에는 A, T, G, C만 사용할 수 있습니다."
            )
            return

        if sequence2:
            sequence2, invalid_bases2 = dna_validation(sequence2)

            if invalid_bases2:
                invalid_text = ", ".join(sorted(invalid_bases2))
                st.error(
                    f"비교 서열에 잘못된 문자가 포함되어 있습니다: "
                    f"{invalid_text}"
                )
                return

        st.session_state.dna_sequence = sequence
        st.session_state.dna_target = target
        st.session_state.dna_sequence2 = sequence2

        go_to("dna_result")
        st.rerun()


# ---------------------------------------------------------
# DNA 결과 화면
# ---------------------------------------------------------

def show_dna_result():
    st.title("🧬 DNA 분석 결과")

    home_button()

    sequence = st.session_state.dna_sequence
    target = st.session_state.dna_target
    sequence2 = st.session_state.dna_sequence2

    counts = count_bases(sequence)
    gc_content = calculate_gc_content(sequence)
    complement = complementary_sequence(sequence)

    st.subheader("기본 분석")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("DNA 길이", f"{len(sequence)} bp")

    with col2:
        st.metric("GC 함량", f"{gc_content:.2f}%")

    with col3:
        st.metric("AT 함량", f"{100 - gc_content:.2f}%")

    st.divider()

    st.subheader("염기 조성")

    base_df = pd.DataFrame({
        "염기": list(counts.keys()),
        "개수": list(counts.values())
    })

    fig = px.bar(
        base_df,
        x="염기",
        y="개수",
        title="DNA 염기별 개수"
    )

    st.plotly_chart(fig, use_container_width=True)

    st.write(
        f"**A:** {counts['A']}개  |  "
        f"**T:** {counts['T']}개  |  "
        f"**G:** {counts['G']}개  |  "
        f"**C:** {counts['C']}개"
    )

    st.divider()

    st.subheader("상보적 DNA 서열")

    st.code(complement)

    st.divider()

    st.subheader("특정 염기서열 검색")

    if target:
        target = target.upper().replace(" ", "")

        invalid_target = set(target) - set("ATGC")

        if invalid_target:
            st.warning(
                "검색할 염기서열에는 A, T, G, C만 입력할 수 있습니다."
            )
        else:
            positions = find_sequence(sequence, target)

            if positions:
                position_text = ", ".join(map(str, positions))

                st.success(
                    f"'{target}' 서열을 {len(positions)}회 찾았습니다."
                )

                st.write(
                    f"시작 위치: {position_text}"
                )

            else:
                st.info(
                    f"'{target}' 서열을 찾을 수 없습니다."
                )

    else:
        st.info("입력 화면에서 검색할 염기서열을 지정할 수 있습니다.")

    st.divider()

    st.subheader("DNA 서열 비교")

    if sequence2:

        differences = compare_sequences(sequence, sequence2)

        if not differences:
            st.success("두 DNA 서열은 동일합니다.")

        else:
            st.warning(
                f"총 {len(differences)}개의 차이가 발견되었습니다."
            )

            difference_df = pd.DataFrame(differences)

            st.dataframe(
                difference_df,
                use_container_width=True,
                hide_index=True
            )

    else:
        st.info("비교할 DNA 서열을 입력하지 않았습니다.")

    st.divider()

    if st.button("새로운 DNA 분석하기", use_container_width=True):
        go_to("dna_input")
        st.rerun()


# ---------------------------------------------------------
# BCA 입력 화면
# ---------------------------------------------------------

def show_bca_input():
    st.title("🧪 BCA 단백질 정량 분석")

    home_button()

    st.write(
        "BCA 실험에서 얻은 표준용액의 농도와 흡광도를 입력하세요."
    )

    st.info(
        "표준용액 데이터를 이용하여 표준곡선을 만들고, "
        "미지 시료의 단백질 농도를 계산합니다."
    )

    st.subheader("표준용액 데이터")

    default_data = pd.DataFrame({
        "단백질 농도 (mg/mL)": [0.0, 0.25, 0.5, 1.0, 1.5, 2.0],
        "흡광도": [0.00, 0.10, 0.21, 0.41, 0.62, 0.82]
    })

    standard_data = st.data_editor(
        default_data,
        num_rows="dynamic",
        use_container_width=True
    )

    st.subheader("미지 시료")

    unknown_absorbance = st.number_input(
        "미지 시료의 흡광도",
        min_value=0.0,
        value=0.30,
        step=0.01
    )

    if st.button(
        "BCA 분석하기",
        type="primary",
        use_container_width=True
    ):

        data = standard_data.dropna()

        if len(data) < 2:
            st.error(
                "표준용액 데이터가 최소 2개 이상 필요합니다."
            )
            return

        concentrations = data["단백질 농도 (mg/mL)"].tolist()
        absorbances = data["흡광도"].tolist()

        if any(value < 0 for value in concentrations):
            st.error("단백질 농도는 음수가 될 수 없습니다.")
            return

        if any(value < 0 for value in absorbances):
            st.error("흡광도는 음수가 될 수 없습니다.")
            return

        if len(set(concentrations)) != len(concentrations):
            st.error(
                "동일한 단백질 농도가 중복되어 있습니다."
            )
            return

        slope, intercept, r_squared = calculate_regression(
            concentrations,
            absorbances
        )

        unknown_concentration = calculate_unknown_concentration(
            unknown_absorbance,
            slope,
            intercept
        )

        if unknown_concentration is None:
            st.error(
                "회귀식을 계산할 수 없습니다. "
                "표준용액 데이터를 확인해주세요."
            )
            return

        st.session_state.bca_data = data
        st.session_state.bca_unknown_absorbance = unknown_absorbance
        st.session_state.bca_slope = slope
        st.session_state.bca_intercept = intercept
        st.session_state.bca_r_squared = r_squared
        st.session_state.bca_unknown_concentration = unknown_concentration

        go_to("bca_result")
        st.rerun()


# ---------------------------------------------------------
# BCA 결과 화면
# ---------------------------------------------------------

def show_bca_result():
    st.title("🧪 BCA 분석 결과")

    home_button()

    data = st.session_state.bca_data
    unknown_absorbance = st.session_state.bca_unknown_absorbance
    slope = st.session_state.bca_slope
    intercept = st.session_state.bca_intercept
    r_squared = st.session_state.bca_r_squared
    unknown_concentration = st.session_state.bca_unknown_concentration

    st.subheader("표준곡선")

    x = data["단백질 농도 (mg/mL)"].to_numpy()
    y = data["흡광도"].to_numpy()

    graph_df = pd.DataFrame({
        "단백질 농도": x,
        "흡광도": y,
        "예측 흡광도": slope * x + intercept
    })

    fig = px.scatter(
        graph_df,
        x="단백질 농도",
        y="흡광도",
        title="BCA 단백질 표준곡선",
        labels={
            "단백질 농도": "단백질 농도 (mg/mL)",
            "흡광도": "흡광도"
        }
    )

    fig.add_scatter(
        x=x,
        y=slope * x + intercept,
        mode="lines",
        name="선형 회귀선"
    )

    st.plotly_chart(fig, use_container_width=True)

    st.divider()

    st.subheader("분석 결과")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric(
            "미지 시료 흡광도",
            f"{unknown_absorbance:.3f}"
        )

    with col2:
        st.metric(
            "단백질 농도",
            f"{unknown_concentration:.3f} mg/mL"
        )

    with col3:
        st.metric(
            "R²",
            f"{r_squared:.4f}"
        )

    st.divider()

    st.subheader("계산에 사용된 회귀식")

    st.code(
        f"y = {slope:.4f}x + {intercept:.4f}"
    )

    st.write(
        "위 회귀식을 이용하여 미지 시료의 단백질 농도를 계산했습니다."
    )

    st.divider()

    st.subheader("표준용액 데이터")

    st.dataframe(
        data,
        use_container_width=True,
        hide_index=True
    )

    st.divider()

    st.subheader("결과 저장")

    result_df = pd.DataFrame({
        "분석 항목": [
            "미지 시료 흡광도",
            "단백질 농도 (mg/mL)",
            "회귀식 기울기",
            "회귀식 절편",
            "R²"
        ],
        "결과": [
            unknown_absorbance,
            unknown_concentration,
            slope,
            intercept,
            r_squared
        ]
    })

    csv_data = result_df.to_csv(
        index=False
    ).encode("utf-8-sig")

    st.download_button(
        label="분석 결과 CSV 저장",
        data=csv_data,
        file_name="BCA_analysis_result.csv",
        mime="text/csv",
        use_container_width=True
    )

    if st.button(
        "새로운 BCA 분석하기",
        use_container_width=True
    ):
        go_to("bca_input")
        st.rerun()


# ---------------------------------------------------------
# 화면 이동
# ---------------------------------------------------------

if st.session_state.page == "home":
    show_home()

elif st.session_state.page == "dna_input":
    show_dna_input()

elif st.session_state.page == "dna_result":
    show_dna_result()

elif st.session_state.page == "bca_input":
    show_bca_input()

elif st.session_state.page == "bca_result":
    show_bca_result()
