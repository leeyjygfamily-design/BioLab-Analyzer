import streamlit as st
import pandas as pd
import plotly.express as px
import requests
from Bio.Align import PairwiseAligner

st.set_page_config(page_title="DNA 분석", page_icon="🧬", layout="wide")

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API = {"tool": "BioLabAnalyzer", "email": st.secrets["NCBI_EMAIL"], "api_key": st.secrets["NCBI_API_KEY"]}

# -------------------- 공통 함수 --------------------
def ncbi_get(endpoint, db, **params):
    r = requests.get(NCBI_URL + endpoint, params={"db": db, **params, **API}, timeout=20)
    r.raise_for_status()
    return r


def search_db(db, term, retmax=10):
    data = ncbi_get("esearch.fcgi", db, term=term, retmode="json", retmax=retmax).json()["esearchresult"]
    return data["idlist"], int(data["count"])


def summaries(db, ids):
    if not ids:
        return []
    data = ncbi_get("esummary.fcgi", db, id=",".join(ids), retmode="json").json()["result"]
    return [data.get(i, {}) for i in ids]


def fetch_sequence(record_id):
    text = ncbi_get("efetch.fcgi", "nuccore", id=record_id, rettype="fasta", retmode="text").text
    return "".join(x.strip() for x in text.splitlines() if not x.startswith(">" )).upper()


def clean_dna(seq):
    return "".join(seq.upper().split())


def base_counts(seq):
    return {b: seq.count(b) for b in "ATGC"}


def gc_percent(seq):
    return 100 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0


# -------------------- 번역 / ORF --------------------
codon_table = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

amino_names = {
    "A":"Alanine","R":"Arginine","N":"Asparagine","D":"Aspartic acid","C":"Cysteine",
    "E":"Glutamic acid","Q":"Glutamine","G":"Glycine","H":"Histidine","I":"Isoleucine",
    "L":"Leucine","K":"Lysine","M":"Methionine","F":"Phenylalanine","P":"Proline",
    "S":"Serine","T":"Threonine","W":"Tryptophan","Y":"Tyrosine","V":"Valine"
}


def translate(seq):
    return "".join(codon_table.get(seq[i:i+3], "?") for i in range(0, len(seq)-2, 3))


def find_orfs(seq):
    orfs, stops = [], {"TAA", "TAG", "TGA"}
    for frame in range(3):
        for start in range(frame, len(seq)-2, 3):
            if seq[start:start+3] != "ATG":
                continue
            for end in range(start+3, len(seq)-2, 3):
                if seq[end:end+3] in stops:
                    dna = seq[start:end+3]
                    orfs.append({"frame": frame+1, "start": start+1, "end": end+3,
                                 "dna": dna, "protein": translate(dna)})
                    break
    return orfs


# -------------------- 서열 정렬 / 변이 분석 --------------------
def align_variants(reference, sample):
    """긴 기준 서열에서 비교 서열과 가장 잘 맞는 구간을 찾아 SNP/삽입/결실을 반환합니다."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1

    alignment = aligner.align(reference, sample)[0]
    coords = alignment.coordinates
    variants, matches, columns = [], 0, 0

    for k in range(len(coords[0]) - 1):
        r1, r2 = int(coords[0][k]), int(coords[0][k+1])
        s1, s2 = int(coords[1][k]), int(coords[1][k+1])
        dr, ds = r2-r1, s2-s1

        if dr and ds:  # 서로 대응되는 구간
            n = min(dr, ds)
            for j in range(n):
                ref, alt = reference[r1+j], sample[s1+j]
                columns += 1
                if ref == alt:
                    matches += 1
                else:
                    variants.append({"기준 위치": r1+j+1, "기준 염기": ref, "비교 염기": alt, "유형": "SNP"})

        elif dr:  # 비교 서열에서 빠진 염기
            columns += dr
            for j in range(dr):
                variants.append({"기준 위치": r1+j+1, "기준 염기": reference[r1+j], "비교 염기": "-", "유형": "Deletion"})

        elif ds:  # 비교 서열에 추가된 염기
            columns += ds
            pos = r1 + 1
            for j in range(ds):
                variants.append({"기준 위치": pos, "기준 염기": "-", "비교 염기": sample[s1+j], "유형": "Insertion"})

    ref_start, ref_end = int(coords[0][0])+1, int(coords[0][-1])
    sample_covered = int(coords[1][-1] - coords[1][0])
    identity = 100 * matches / columns if columns else 0
    coverage = 100 * sample_covered / len(sample) if sample else 0

    return {
        "variants": variants,
        "start": ref_start,
        "end": ref_end,
        "identity": identity,
        "coverage": coverage,
        "score": alignment.score
    }


def substitution_type(ref, alt):
    if ref == "-" or alt == "-":
        return "-"
    return "Transition" if (ref, alt) in {("A","G"),("G","A"),("C","T"),("T","C")} else "Transversion"


def add_orf_effect(variant, reference, selected_orf):
    """SNP가 선택한 ORF 안에 있을 때만 코돈/아미노산 변화를 계산합니다."""
    if variant["유형"] != "SNP" or not selected_orf:
        return variant

    pos = variant["기준 위치"]
    if not (selected_orf["start"] <= pos <= selected_orf["end"]):
        return variant

    offset = pos - selected_orf["start"]
    codon_start = selected_orf["start"] - 1 + (offset // 3) * 3
    ref_codon = reference[codon_start:codon_start+3]
    alt_codon = list(ref_codon)
    alt_codon[offset % 3] = variant["비교 염기"]
    alt_codon = "".join(alt_codon)
    ref_aa, alt_aa = codon_table.get(ref_codon, "?"), codon_table.get(alt_codon, "?")

    if ref_aa == alt_aa:
        effect = "Synonymous"
    elif alt_aa == "*":
        effect = "Nonsense"
    else:
        effect = "Missense"

    return {**variant, "기준 코돈": ref_codon, "변이 코돈": alt_codon,
            "기준 아미노산": ref_aa, "변이 아미노산": alt_aa, "영향": effect}


def save_alignment(sample):
    result = align_variants(st.session_state.sequence, sample)
    st.session_state.variant_sample = sample
    st.session_state.variant_result = result


# ==================== 화면 ====================
st.title("🧬 DNA Analysis")
st.write("NCBI 유전자 서열을 가져와 DNA 특성, ORF, 아미노산 및 변이를 분석합니다.")

# ① NCBI 검색
st.subheader("① NCBI 유전자 검색")
c1, c2 = st.columns(2)
gene = c1.text_input("유전자명", placeholder="예: TP53")
organism = c2.text_input("생물 종", value="Homo sapiens")

if st.button("NCBI에서 검색하기", type="primary", use_container_width=True):
    if not gene:
        st.warning("검색할 유전자명을 입력해주세요.")
    else:
        try:
            ids, total = search_db("nuccore", f"{gene}[Gene Name] AND {organism}[Organism]", 5)
            rows = summaries("nuccore", ids)
            st.session_state.ncbi_total = total
            st.session_state.ncbi_results = [
                {"ID": i, "제목": r.get("title", "정보 없음"), "길이": r.get("slen", "정보 없음")}
                for i, r in zip(ids, rows)
            ]
        except requests.RequestException:
            st.error("NCBI 서버와 연결할 수 없습니다.")

if "ncbi_results" in st.session_state:
    results = st.session_state.ncbi_results
    st.write(f"검색 결과: 총 {st.session_state.get('ncbi_total', len(results)):,}개 (상위 {len(results)}개 표시)")
    st.dataframe(pd.DataFrame(results), use_container_width=True, hide_index=True)
    selected = st.selectbox("분석할 서열", results, format_func=lambda x: f"{x['ID']} | {x['제목']}")

    if st.button("이 서열 분석하기", use_container_width=True):
        try:
            st.session_state.sequence = fetch_sequence(selected["ID"])
            for key in ("orfs", "variant_result", "variant_sample"):
                st.session_state.pop(key, None)
            st.success("DNA 서열을 가져왔습니다.")
        except requests.RequestException:
            st.error("DNA 서열을 가져오는 데 실패했습니다.")


if "sequence" not in st.session_state:
    st.info("먼저 위에서 기준 DNA 서열을 선택해주세요.")
    st.stop()

sequence = st.session_state.sequence

# ② 기본 분석
st.divider()
st.subheader("② DNA 기본 분석")
c1, c2 = st.columns(2)
c1.metric("DNA 길이", f"{len(sequence):,} bp")
c2.metric("GC 함량", f"{gc_percent(sequence):.2f}%")
counts = base_counts(sequence)
fig = px.bar(pd.DataFrame({"염기": counts.keys(), "개수": counts.values()}), x="염기", y="개수", title="DNA 염기 조성")
st.plotly_chart(fig, use_container_width=True)
with st.expander("GC 함량은 무엇을 의미하나요?"):
    st.write("GC 함량은 전체 염기 중 G와 C가 차지하는 비율입니다. 서열의 특성을 비교하는 기본 지표로 활용됩니다.")

# ③ ORF
st.divider()
st.subheader("③ ORF 및 아미노산 분석")
with st.expander("ORF란?"):
    st.write("ORF(Open Reading Frame)는 시작 코돈에서 종결 코돈까지 이어져 단백질로 번역될 가능성이 있는 염기서열 구간입니다.")

if st.button("ORF 분석하기", type="primary", use_container_width=True):
    st.session_state.orfs = find_orfs(sequence)

selected_orf = None
if "orfs" in st.session_state:
    orfs = st.session_state.orfs
    if not orfs:
        st.warning("완전한 ORF를 찾지 못했습니다.")
    else:
        st.success(f"{len(orfs)}개의 ORF를 찾았습니다.")
        labels = [f"ORF {i+1} | Frame {o['frame']} | {o['start']}–{o['end']}" for i, o in enumerate(orfs)]
        idx = st.selectbox("분석할 ORF", range(len(orfs)), format_func=lambda i: labels[i])
        selected_orf = orfs[idx]
        st.session_state.selected_orf = selected_orf

        c1, c2 = st.columns(2)
        c1.metric("ORF 길이", f"{len(selected_orf['dna'])} bp")
        c2.metric("아미노산 길이", f"{max(0, len(selected_orf['protein'])-1)} aa")
        st.write("번역된 아미노산 서열")
        st.code(selected_orf["protein"])

        options = {f"{name} ({code})": code for code, name in amino_names.items()}
        amino_label = st.selectbox("검색할 아미노산", options)
        code = options[amino_label]
        positions = [i+1 for i, aa in enumerate(selected_orf["protein"]) if aa == code]
        st.write(f"{amino_label}: {len(positions)}개")
        if positions:
            st.write("위치:", ", ".join(map(str, positions)))
else:
    selected_orf = st.session_state.get("selected_orf")

# ④ 변이 분석
st.divider()
st.subheader("④ DNA 변이 분석")
st.write("비교 서열이 더 짧거나 길어도 기준 DNA에서 가장 유사한 구간을 찾아 정렬한 뒤 변이를 분석합니다.")

source = st.radio("비교 데이터", ["직접 DNA 입력", "NCBI 비교 서열", "ClinVar 알려진 변이"], horizontal=True)

if source == "직접 DNA 입력":
    sample_text = st.text_area("비교할 DNA 서열", height=130, placeholder="ATGCGT...")
    if st.button("변이 분석하기", type="primary", use_container_width=True, key="manual_align"):
        sample = clean_dna(sample_text)
        if not sample:
            st.warning("비교할 DNA를 입력해주세요.")
        elif any(b not in "ATGC" for b in sample):
            st.error("A, T, G, C만 입력해주세요.")
        else:
            save_alignment(sample)

elif source == "NCBI 비교 서열":
    query = st.text_input("유전자명, accession 또는 검색어", placeholder="예: TP53 Homo sapiens")
    if st.button("NCBI 비교 서열 검색", use_container_width=True):
        try:
            ids, total = search_db("nuccore", query, 10)
            rows = summaries("nuccore", ids)
            st.session_state.compare_total = total
            st.session_state.compare_results = [
                {"ID": i, "제목": r.get("title", "제목 없음"), "Accession": r.get("accessionversion", ""), "길이": r.get("slen", 0)}
                for i, r in zip(ids, rows)
            ]
        except requests.RequestException:
            st.error("NCBI 검색에 실패했습니다.")

    if "compare_results" in st.session_state:
        rows = st.session_state.compare_results
        st.write(f"검색 결과: 총 {st.session_state.get('compare_total', len(rows)):,}개 (상위 {len(rows)}개 표시)")
        choice = st.selectbox("비교할 서열", rows, format_func=lambda x: f"{x['Accession']} | {x['제목']} | {x['길이']:,} bp")
        if st.button("선택한 서열 변이 분석하기", type="primary", use_container_width=True):
            try:
                save_alignment(fetch_sequence(choice["ID"]))
            except requests.RequestException:
                st.error("비교 서열을 가져오지 못했습니다.")

else:  # ClinVar
    query = st.text_input("유전자명 또는 변이", placeholder="예: BRCA1 또는 TP53")
    if st.button("ClinVar 검색", use_container_width=True):
        try:
            ids, total = search_db("clinvar", query, 10)
            rows = summaries("clinvar", ids)
            st.session_state.clinvar_total = total
            st.session_state.clinvar_results = [{"ID": i, **r} for i, r in zip(ids, rows)]
        except requests.RequestException:
            st.error("ClinVar 검색에 실패했습니다.")

    if "clinvar_results" in st.session_state:
        rows = st.session_state.clinvar_results
        st.write(f"검색 결과: 총 {st.session_state.get('clinvar_total', len(rows)):,}개 (상위 {len(rows)}개 표시)")
        choice = st.selectbox("분석할 ClinVar 변이", rows, format_func=lambda x: f"{x['ID']} | {x.get('title', '정보 없음')}")

        if st.button("선택한 변이 분석하기", type="primary", use_container_width=True):
            st.session_state.selected_clinvar = choice

    if "selected_clinvar" in st.session_state:
        v = st.session_state.selected_clinvar
        classification = v.get("germline_classification", {})
        if not isinstance(classification, dict):
            classification = {}
        genes = v.get("genes", [])
        gene_text = ", ".join(g.get("symbol", "") for g in genes if isinstance(g, dict)) if isinstance(genes, list) else "정보 없음"

        st.subheader("ClinVar 변이 분석 결과")
        st.write("**변이:**", v.get("title", "정보 없음"))
        st.write("**유전자:**", gene_text or "정보 없음")
        st.write("**변이 유형:**", v.get("variant_type", v.get("obj_type", "정보 없음")))
        st.write("**임상적 분류:**", classification.get("description", "정보 없음"))
        st.caption("ClinVar는 알려진 변이 기록을 조회하는 기능이며, 전체 비교 DNA 서열을 제공하는 기능과는 다릅니다.")

# 정렬 결과는 직접 입력/NCBI 비교에서 표시
if source != "ClinVar 알려진 변이" and "variant_result" in st.session_state:
    r = st.session_state.variant_result
    variants = r["variants"]

    st.divider()
    st.subheader("변이 분석 결과")
    c1, c2, c3 = st.columns(3)
    c1.metric("가장 유사한 기준 구간", f"{r['start']:,}–{r['end']:,} bp")
    c2.metric("서열 일치도", f"{r['identity']:.2f}%")
    c3.metric("비교 서열 정렬 범위", f"{r['coverage']:.2f}%")

    if not variants:
        st.success("정렬된 구간에서 변이가 발견되지 않았습니다.")
    else:
        selected_orf = st.session_state.get("selected_orf")
        detailed = []
        for v in variants:
            item = {**v, "치환 유형": substitution_type(v["기준 염기"], v["비교 염기"])}
            detailed.append(add_orf_effect(item, sequence, selected_orf))

        df = pd.DataFrame(detailed)
        st.write(f"총 {len(df)}개의 차이를 찾았습니다.")
        st.dataframe(df, use_container_width=True, hide_index=True)

        counts = df["유형"].value_counts().rename_axis("유형").reset_index(name="개수")
        st.plotly_chart(px.bar(counts, x="유형", y="개수", title="변이 유형"), use_container_width=True)

        positions = pd.DataFrame({"위치": [v["기준 위치"] for v in variants], "값": 1, "유형": [v["유형"] for v in variants]})
        fig = px.scatter(positions, x="위치", y="값", hover_name="유형", title="기준 DNA 내 변이 위치")
        fig.update_yaxes(visible=False)
        st.plotly_chart(fig, use_container_width=True)
