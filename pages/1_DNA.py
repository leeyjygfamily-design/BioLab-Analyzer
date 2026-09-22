import streamlit as st
import pandas as pd
import plotly.express as px


st.set_page_config(
    page_title="DNA 분석",
    page_icon="🧬"
)


# DNA 서열을 검사하는 함수
def validate_dna(sequence):
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    valid_bases = set("ATGC")
    invalid_bases = set(sequence) - valid_bases

    return sequence, invalid_bases


# DNA의 염기 개수를 세는 함수
def count_bases(sequence):
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }


# GC 함량을 계산하는 함수
def calculate_gc(sequence):
    if len(sequence) == 0:
        return 0

    gc = sequence.count("G") + sequence.count("C")

    return gc / len(sequence) * 100


# 상보적 염기서열을 만드는 함수
def make_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    return "".join(complement[base] for base in sequence)


# 특정 염기서열의 위치를 찾는 함수
def find_sequence(sequence, target):
    positions = []
    start = 0

    while True:
        position = sequence.find(target, start)

        if position == -1:
            break

        positions.append(position + 1)
        start = position + 1

    return positions


# --------------------------------------------------
# 화면
# --------------------------------------------------

st.title("🧬 DNA 서열 분석")

st.write(
    "DNA 서열을 입력하면 염기 조성, GC 함량, "
    "상보적 서열 등을 분석합니다."
)

sequence = st.text_area(
    "DNA 서열 입력",
    placeholder="예: ATGCGTACCGTA"
)

target = st.text_input(
    "검색할 염기서열",
    placeholder="예: ATG"
)

if st.button("DNA 분석하기", type="primary"):

    sequence, invalid_bases = validate_dna(sequence)

    if not sequence:
        st.error("DNA 서열을 입력해주세요.")

    elif invalid_bases:
        st.error(
            "A, T, G, C 이외의 문자가 포함되어 있습니다."
        )

    else:

        counts = count_bases(sequence)
        gc_content = calculate_gc(sequence)
        complement = make_complement(sequence)

        st.divider()

        st.subheader("분석 결과")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                "DNA 길이",
                f"{len(sequence)} bp"
            )

        with col2:
            st.metric(
                "GC 함량",
                f"{gc_content:.2f}%"
            )

        with col3:
            st.metric(
                "AT 함량",
                f"{100 - gc_content:.2f}%"
            )

        st.subheader("염기 조성")

        base_data = pd.DataFrame({
            "염기": list(counts.keys()),
            "개수": list(counts.values())
        })

        fig = px.bar(
            base_data,
            x="염기",
            y="개수",
            title="DNA 염기별 개수"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

        st.subheader("상보적 DNA 서열")

        st.code(complement)

        st.subheader("특정 염기서열 검색")

        if target:

            target = target.upper().replace(" ", "")

            if set(target) - set("ATGC"):
                st.warning(
                    "검색어에는 A, T, G, C만 사용할 수 있습니다."
                )

            else:

                positions = find_sequence(
                    sequence,
                    target
                )

                if positions:
                    st.success(
                        f"'{target}' 서열을 "
                        f"{len(positions)}회 찾았습니다."
                    )

                    st.write(
                        "시작 위치:",
                        ", ".join(map(str, positions))
                    )

                else:
                    st.info(
                        f"'{target}' 서열을 찾지 못했습니다."
                    )
