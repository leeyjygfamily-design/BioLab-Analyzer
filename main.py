import streamlit as st

# 페이지 설정
st.set_page_config(
    page_title="BioLab Analyzer",
    page_icon="🧬",
    layout="wide"
)

# DNA 분석 페이지
dna_page = st.Page(
    "pages/1_DNA.py",
    title="DNA 분석",
    icon="🧬"
)

# BCA 단백질 정량 페이지
bca_page = st.Page(
    "pages/2_BCA_protein.py",
    title="BCA 단백질 정량",
    icon="🧪"
)

# PCR 프라이머 분석 페이지
pcr_page = st.Page(
    "pages/3_PCR.py",
    title="PCR 프라이머 분석",
    icon="🔬"
)

# AI 챗봇 페이지
ai_page = st.Page(
    "pages/4_AI_teacher.py",
    title="AI 생명공학 연구 도우미",
    icon="🤖"
)

# 사이드바 페이지 이동 메뉴
pg = st.navigation(
    {
        "BioLab Analyzer": [
            dna_page,
            bca_page,
            pcr_page,
            ai_page
        ]
    }
)

# 현재 선택된 페이지 실행
pg.run()
