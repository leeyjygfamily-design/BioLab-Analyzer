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
# ORF 및 아미노산 분석
# --------------------------------------------------

# DNA 코돈을 아미노산으로 변환하기 위한 표입니다.
codon_table = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",

    "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*",

    "TGT": "C", "TGC": "C",
    "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L",

    "CCT": "P", "CCC": "P",
    "CCA": "P", "CCG": "P",

    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",

    "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I",
    "ATA": "I", "ATG": "M",

    "ACT": "T", "ACC": "T",
    "ACA": "T", "ACG": "T",

    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",

    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V",

    "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A",

    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",

    "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G"
}


# 아미노산 한 글자를 이름으로 바꾸기 위한 표입니다.
amino_acid_names = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine"
}


def find_orfs(sequence):
    """
    DNA 서열의 3가지 reading frame에서
    시작 코돈 ATG부터 종결 코돈까지 ORF를 찾습니다.
    """

    orfs = []

    stop_codons = {"TAA", "TAG", "TGA"}

    # DNA를 읽기 시작하는 위치를 0, 1, 2로 바꿔봅니다.
    for frame in range(3):

        i = frame

        while i <= len(sequence) - 3:

            codon = sequence[i:i + 3]

            # 시작 코돈을 발견하면 ORF 탐색을 시작합니다.
            if codon == "ATG":

                start = i
                j = i + 3

                while j <= len(sequence) - 3:

                    current_codon = sequence[j:j + 3]

                    # 종결 코돈을 찾으면 ORF를 완성합니다.
                    if current_codon in stop_codons:

                        end = j + 3

                        dna_sequence = sequence[start:end]

                        amino_sequence = translate_dna(
                            dna_sequence
                        )

                        orfs.append({
                            "reading_frame": frame + 1,
                            "start": start + 1,
                            "end": end,
                            "length": end - start,
                            "dna": dna_sequence,
                            "protein": amino_sequence
                        })

                        break

                    j += 3

                # 다음 ATG를 찾기 위해 계속 이동합니다.

            i += 3

    return orfs


def translate_dna(sequence):
    """
    DNA 서열을 3개씩 읽어서 아미노산 서열로 변환합니다.
    """

    protein = []

    for i in range(0, len(sequence) - 2, 3):

        codon = sequence[i:i + 3]

        amino_acid = codon_table.get(
            codon,
            "?"
        )

        protein.append(amino_acid)

    return "".join(protein)

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
# --------------------------------------------------
# ORF 분석
# --------------------------------------------------

st.divider()

st.subheader("🧬 ORF 및 아미노산 분석")

st.write(
    "DNA 서열에서 시작 코돈(ATG)부터 "
    "종결 코돈(TAA, TAG, TGA)까지의 "
    "ORF를 탐색하고 아미노산 서열로 변환합니다."
)

if st.button(
    "ORF 찾기",
    use_container_width=True
):

    orfs = find_orfs(sequence)

    if not orfs:

        st.warning(
            "입력된 DNA 서열에서 "
            "완전한 ORF를 찾지 못했습니다."
        )

    else:

        st.success(
            f"{len(orfs)}개의 ORF를 찾았습니다."
        )

        orf_data = []

        for index, orf in enumerate(orfs):

            orf_data.append({
                "ORF": f"ORF {index + 1}",
                "Reading frame": orf["reading_frame"],
                "시작 위치": orf["start"],
                "종료 위치": orf["end"],
                "DNA 길이": orf["length"],
                "아미노산 길이": len(
                    orf["protein"]
                ) - 1
            })

        orf_df = pd.DataFrame(
            orf_data
        )

        st.dataframe(
            orf_df,
            use_container_width=True,
            hide_index=True
        )

        # ORF 하나를 선택합니다.
        selected_orf_number = st.selectbox(
            "자세히 볼 ORF를 선택하세요.",
            range(
                1,
                len(orfs) + 1
            )
        )

        selected_orf = orfs[
            selected_orf_number - 1
        ]

        st.subheader(
            f"ORF {selected_orf_number} 상세 정보"
        )

        st.write(
            f"Reading frame: "
            f"{selected_orf['reading_frame']}"
        )

        st.write(
            f"DNA 위치: "
            f"{selected_orf['start']} ~ "
            f"{selected_orf['end']}"
        )

        st.write(
            f"DNA 길이: "
            f"{selected_orf['length']} bp"
        )

        st.write("DNA 서열")

        st.code(
            selected_orf["dna"]
        )

        st.write("번역된 아미노산 서열")

        st.code(
            selected_orf["protein"]
        )


        # --------------------------------------------------
        # 아미노산 검색
        # --------------------------------------------------

        st.subheader(
            "🔎 특정 아미노산 검색"
        )

        amino_options = {
            f"{name} ({code})": code
            for code, name
            in amino_acid_names.items()
        }

        selected_amino_name = st.selectbox(
            "검색할 아미노산을 선택하세요.",
            list(amino_options.keys())
        )

        selected_amino = amino_options[
            selected_amino_name
        ]

        protein_sequence = (
            selected_orf["protein"]
        )

        amino_positions = []

        for index, amino in enumerate(
            protein_sequence
        ):

            # 종결 코돈은 아미노산 검색에서 제외합니다.
            if amino == selected_amino:

                amino_positions.append(
                    index + 1
                )

        if amino_positions:

            st.success(
                f"{selected_amino_name}이 "
                f"{len(amino_positions)}개 발견되었습니다."
            )

            st.write(
                "아미노산 위치:",
                ", ".join(
                    map(
                        str,
                        amino_positions
                    )
                )
            )

            # 해당 아미노산의 DNA 코돈도 보여줍니다.
            codon_results = []

            for position in amino_positions:

                dna_index = (
                    (position - 1) * 3
                )

                codon = selected_orf[
                    "dna"
                ][
                    dna_index:dna_index + 3
                ]

                codon_results.append({
                    "아미노산 위치": position,
                    "아미노산": selected_amino,
                    "DNA 코돈": codon,
                    "DNA 위치": (
                        selected_orf["start"]
                        + dna_index
                    )
                })

            amino_df = pd.DataFrame(
                codon_results
            )

            st.dataframe(
                amino_df,
                use_container_width=True,
                hide_index=True
            )

        else:

            st.info(
                f"선택한 ORF에서 "
                f"{selected_amino_name}을 "
                "찾지 못했습니다."
            )
