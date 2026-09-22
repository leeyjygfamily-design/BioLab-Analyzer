import streamlit as st

# 페이지 설정
st.set_page_config(
    page_title="BioLab Analyzer",
    page_icon="🧬",
    layout="wide"
)

# 분석 페이지 등록
dna_page = st.Page(
    "pages/1_DNA.py",
    title="DNA 분석",
    icon="🧬"
)

bca_page = st.Page(
    "pages/2_BCA_protein.py",
    title="BCA 단백질 정량",
    icon="🧪"
)

# 페이지 이동 메뉴 만들기
pg = st.navigation(
    {
        "BioLab Analyzer": [
            dna_page,
            bca_page
        ]
    }
)

# 현재 선택된 페이지 실행
pg.run()
