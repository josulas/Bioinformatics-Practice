import xml.etree.ElementTree as ET


FILE_NAME = "blast_output_remote_1.xml"


def get_list_match(root, output: list, tag_name: str):
    if len(root):
        finding = root.findall(tag_name)
        if len(finding):
            output.extend(finding)
        else:
            for child in root:
                get_list_match(child, output, tag_name)


def get_dict_tags_and_seqs(root):
    hits = []
    get_list_match(root, hits, "Hit")
    list_output = []
    for hit in hits:
        id = hit.find('Hit_def').text
        seqs = []
        scores = []
        get_list_match(root, seqs, 'Hsp_hseq')
        get_list_match(root, scores, 'Hsp_score')
        for index, seq in enumerate(seqs):
            list_output.append((id + f"; hsp_index: {index}", seq.text.replace('-', ''), int(scores[index].text)))
    return list_output

tree = ET.parse(FILE_NAME)
root = tree.getroot()
sequences = get_dict_tags_and_seqs(root)
sequences.sort(key= lambda x: x[2], reverse=True)
with open('best10.fasta', 'w', encoding='UTF-8') as file:
    for i in range(10):
        file.write(f"> {sequences[i][0]}\n")
        file.write(f"{sequences[i][1]}\n\n")
