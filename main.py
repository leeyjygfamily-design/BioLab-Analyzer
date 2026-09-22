import streamlit as st

st.set_page_config(
    page_title="BioLab Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 BioLab Analyzer")
st.subheader("생체분자 데이터 분석 프로그램")

st.write(
    "생명공학 실험에서 얻은 DNA 서열과 "
    "BCA 단백질 정량 데이터를 분석할 수 있습니다."
)

st.divider()

col1, col2 = st.columns(2)

with col1:
    st.subheader("🧬 DNA 서열 분석")
    st.write(
        "DNA 염기 조성, GC 함량, "
        "상보적 서열 등을 분석합니다."
    )

    if st.button("DNA 분석 시작 →", use_container_width=True):
        st.switch_page("pages/1_DNA.py")


with col2:
    st.subheader("🧪 BCA 단백질 정량")
    st.write(
        "BCA 실험 데이터를 이용하여 "
        "미지 시료의 단백질 농도를 계산합니다."
    )

    if st.button("BCA 분석 시작 →", use_container_width=True):
        st.switch_page("pages/2_BCA_protein.py")
