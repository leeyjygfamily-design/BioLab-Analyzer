import streamlit as st
import pandas as pd
import plotly.express as px
import requests


# --------------------------------------------------
# 기본 설정
# --------------------------------------------------

st.set_page_config(
    page_title="DNA 분석",
    page_icon="🧬",
    layout="wide"
)


# --------------------------------------------------
# NCBI API 설정
# --------------------------------------------------

NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Streamlit Secrets에서 API 키와 이메일을 가져옵니다.
try:
    NCBI_API_KEY = st.secrets["NCBI_API_KEY"]
    NCBI_EMAIL = st.secrets["NCBI_EMAIL"]
    NCBI_READY = True

except Exception:
    NCBI_API_KEY = ""
    NCBI_EMAIL = ""
    NCBI_READY = False


# --------------------------------------------------
# DNA 분석 함수
# --------------------------------------------------

def validate_dna(sequence):
    """DNA 서열에 잘못된 문자가 있는지 확인합니다."""

    sequence = (
        sequence
        .upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
    )

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


def calculate_gc(sequence):
    """DNA 서열의 GC 함량을 계산합니다."""

    if len(sequence) == 0:
        return 0

    gc_count = (
        sequence.count("G")
        + sequence.count("C")
    )

    return gc_count / len(sequence) * 100


def make_complement(sequence):
    """DNA의 상보적 염기서열을 만듭니다."""

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    return "".join(
        complement[base]
        for base in sequence
    )


def find_sequence(sequence, target):
    """특정 염기서열이 나타나는 위치를 찾습니다."""

    positions = []

    start = 0

    while True:

        position = sequence.find(
            target,
            start
        )

        if position == -1:
            break

        # 위치는 사람이 보기 쉽게 1부터 표시합니다.
        positions.append(position + 1)

        start = position + 1

    return positions


def compare_sequences(sequence1, sequence2):
    """두 DNA 서열의 차이를 비교합니다."""

    max_length = max(
        len(sequence1),
        len(sequence2)
    )

    differences = []

    for i in range(max_length):

        base1 = (
            sequence1[i]
            if i < len(sequence1)
            else "-"
        )

        base2 = (
            sequence2[i]
            if i < len(sequence2)
            else "-"
        )

        if base1 != base2:

            differences.append({
                "위치": i + 1,
                "서열 1": base1,
                "서열 2": base2
            })

    return differences


# --------------------------------------------------
# NCBI 검색 함수
# --------------------------------------------------

def search_ncbi(query, organism, max_results=10):
    """NCBI nucleotide 데이터베이스에서 DNA 서열을 검색합니다."""

    search_term = (
        f"{query}[Gene Name] "
        f"AND {organism}[Organism]"
    )

    params = {
        "db": "nuccore",
        "term": search_term,
        "retmode": "json",
        "retmax": max_results,
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_BASE_URL + "esearch.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    data = response.json()

    id_list = data.get(
        "esearchresult",
        {}
    ).get(
        "idlist",
        []
    )

    return id_list


def get_ncbi_summary(id_list):
    """검색된 DNA 기록의 기본 정보를 가져옵니다."""

    if not id_list:
        return []

    params = {
        "db": "nuccore",
        "id": ",".join(id_list),
        "retmode": "json",
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_BASE_URL + "esummary.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    data = response.json()

    result = data.get(
        "result",
        {}
    )

    summaries = []

    for record_id in id_list:

        record = result.get(
            record_id,
            {}
        )

        summaries.append({
            "ID": record_id,
            "제목": record.get(
                "title",
                "정보 없음"
            ),
            "길이": record.get(
                "slen",
                "정보 없음"
            )
        })

    return summaries


def get_ncbi_sequence(record_id):
    """선택한 NCBI 기록의 DNA 서열을 FASTA 형식으로 가져옵니다."""

    params = {
        "db": "nuccore",
        "id": record_id,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_BASE_URL + "efetch.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    fasta_text = response.text.strip()

    # FASTA의 첫 번째 줄은 서열 설명이므로 제거합니다.
    fasta_lines = fasta_text.splitlines()

    sequence_lines = [
        line.strip()
        for line in fasta_lines
        if not line.startswith(">")
    ]

    sequence = "".join(sequence_lines)

    return sequence, fasta_text


# --------------------------------------------------
# NCBI 검색 화면
# --------------------------------------------------

st.title("🧬 DNA 서열 분석")

st.write(
    "DNA 서열을 직접 입력하거나 "
    "NCBI에서 실제 유전자 서열을 검색하여 분석할 수 있습니다."
)


st.divider()

st.subheader("🔎 NCBI에서 유전자 검색")

if not NCBI_READY:

    st.error(
        "NCBI API 설정을 찾을 수 없습니다."
    )

    st.info(
        "Streamlit Secrets에 "
        "NCBI_API_KEY와 NCBI_EMAIL을 등록해주세요."
    )

else:

    col1, col2 = st.columns(2)

    with col1:

        gene_name = st.text_input(
            "유전자명",
            placeholder="예: TP53"
        )

    with col2:

        organism = st.text_input(
            "생물 종",
            value="Homo sapiens",
            placeholder="예: Homo sapiens"
        )

    if st.button(
        "NCBI에서 검색하기",
        type="primary",
        use_container_width=True
    ):

        if not gene_name.strip():

            st.warning(
                "검색할 유전자명을 입력해주세요."
            )

        elif not organism.strip():

            st.warning(
                "생물 종을 입력해주세요."
            )

        else:

            with st.spinner(
                "NCBI에서 DNA 서열을 검색하고 있습니다..."
            ):

                try:

                    id_list = search_ncbi(
                        gene_name.strip(),
                        organism.strip()
                    )

                    if not id_list:

                        st.warning(
                            "검색 결과가 없습니다. "
                            "유전자명과 생물 종을 확인해주세요."
                        )

                    else:

                        summaries = get_ncbi_summary(
                            id_list
                        )

                        st.session_state.ncbi_results = summaries

                        st.success(
                            f"{len(summaries)}개의 검색 결과를 찾았습니다."
                        )

                except requests.exceptions.Timeout:

                    st.error(
                        "NCBI 서버의 응답 시간이 초과되었습니다."
                    )

                except requests.exceptions.RequestException as e:

                    st.error(
                        f"NCBI API 요청 중 오류가 발생했습니다: {e}"
                    )

                except Exception as e:

                    st.error(
                        f"검색 중 오류가 발생했습니다: {e}"
                    )


# --------------------------------------------------
# NCBI 검색 결과
# --------------------------------------------------

if "ncbi_results" in st.session_state:

    st.divider()

    st.subheader("📋 NCBI 검색 결과")

    results = st.session_state.ncbi_results

    result_df = pd.DataFrame(results)

    st.dataframe(
        result_df,
        use_container_width=True,
        hide_index=True
    )

    result_options = {
        f"{item['ID']} | {item['제목']}": item["ID"]
        for item in results
    }

    selected_label = st.selectbox(
        "분석할 DNA 기록을 선택하세요.",
        list(result_options.keys())
    )

    selected_id = result_options[
        selected_label
    ]

    if st.button(
        "선택한 DNA 서열 가져오기",
        use_container_width=True
    ):

        with st.spinner(
            "NCBI에서 DNA 서열을 가져오는 중입니다..."
        ):

            try:

                sequence, fasta_text = (
                    get_ncbi_sequence(
                        selected_id
                    )
                )

                sequence, invalid_bases = (
                    validate_dna(sequence)
                )

                if invalid_bases:

                    st.error(
                        "NCBI에서 가져온 서열에 "
                        "예상하지 못한 문자가 포함되어 있습니다."
                    )

                else:

                    st.session_state.ncbi_sequence = sequence
                    st.session_state.ncbi_fasta = fasta_text
                    st.session_state.ncbi_id = selected_id

                    st.success(
                        f"DNA 서열을 가져왔습니다. "
                        f"총 {len(sequence):,} bp입니다."
                    )

            except requests.exceptions.Timeout:

                st.error(
                    "DNA 서열을 가져오는 시간이 초과되었습니다."
                )

            except requests.exceptions.RequestException as e:

                st.error(
                    f"DNA 서열을 가져오는 중 오류가 발생했습니다: {e}"
                )


# --------------------------------------------------
# NCBI에서 가져온 서열
# --------------------------------------------------

if "ncbi_sequence" in st.session_state:

    st.divider()

    st.subheader("📥 NCBI에서 가져온 DNA 서열")

    st.write(
        f"NCBI ID: "
        f"`{st.session_state.ncbi_id}`"
    )

    st.text_area(
        "가져온 서열",
        value=st.session_state.ncbi_sequence,
        height=150
    )

    with st.expander(
        "FASTA 원본 확인"
    ):

        st.code(
            st.session_state.ncbi_fasta,
            language="text"
        )

    if st.button(
        "이 서열을 분석하기",
        type="primary",
        use_container_width=True
    ):

        st.session_state.dna_sequence = (
            st.session_state.ncbi_sequence
        )

        st.success(
            "NCBI 서열을 아래 분석에 적용했습니다."
        )


# --------------------------------------------------
# 직접 DNA 서열 입력
# --------------------------------------------------

st.divider()

st.subheader("✏️ DNA 서열 직접 입력")

sequence = st.text_area(
    "DNA 서열 입력",
    value=st.session_state.get(
        "dna_sequence",
        ""
    ),
    placeholder="예: ATGCGTACCGTA",
    height=150
)

target = st.text_input(
    "검색할 염기서열",
    placeholder="예: ATG"
)


# --------------------------------------------------
# DNA 분석
# --------------------------------------------------

if st.button(
    "DNA 분석하기",
    type="primary",
    use_container_width=True
):

    sequence, invalid_bases = (
        validate_dna(sequence)
    )

    if not sequence:

        st.error(
            "DNA 서열을 입력해주세요."
        )

    elif invalid_bases:

        invalid_text = ", ".join(
            sorted(invalid_bases)
        )

        st.error(
            f"잘못된 문자가 포함되어 있습니다: "
            f"{invalid_text}"
        )

        st.info(
            "DNA 서열에는 A, T, G, C만 사용할 수 있습니다."
        )

    else:

        counts = count_bases(sequence)

        gc_content = calculate_gc(
            sequence
        )

        complement = make_complement(
            sequence
        )

        st.divider()

        st.subheader("📊 DNA 분석 결과")

        col1, col2, col3 = st.columns(3)

        with col1:

            st.metric(
                "DNA 길이",
                f"{len(sequence):,} bp"
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

        st.divider()

        st.subheader("염기 조성")

        base_data = pd.DataFrame({
            "염기": list(
                counts.keys()
            ),
            "개수": list(
                counts.values()
            )
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

        st.write(
            f"**A:** {counts['A']:,}개  |  "
            f"**T:** {counts['T']:,}개  |  "
            f"**G:** {counts['G']:,}개  |  "
            f"**C:** {counts['C']:,}개"
        )

        st.divider()

        st.subheader("상보적 DNA 서열")

        st.code(
            complement
        )

        st.divider()

        st.subheader("특정 염기서열 검색")

        if target:

            target = (
                target
                .upper()
                .replace(" ", "")
            )

            invalid_target = (
                set(target) - set("ATGC")
            )

            if invalid_target:

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
                        "시작 위치: "
                        + ", ".join(
                            map(
                                str,
                                positions
                            )
                        )
                    )

                else:

                    st.info(
                        f"'{target}' 서열을 "
                        "찾지 못했습니다."
                    )

        else:

            st.info(
                "검색할 염기서열을 입력하면 "
                "해당 서열의 위치를 찾을 수 있습니다."
            )
