# This script contains function to parse parent pathway for Gene Ontology pathway term.
# Pathway term must be in the form of lower cases, without any exclamation marks.

from goatools.obo_parser import GODag
obofile = './go.obo'  # Replace with the path to your Gene Ontology file
go =GODag(obofile)


def parse_parents(go, term_name, og_id=None, parent_path_lst=None):
    if parent_path_lst is None:
        parent_path_lst = []
    for term in go.values():
        if term.name.lower().replace('-',' ') == term_name.lower():
            pathway_term = term
            og_id = term.id
    if pathway_term is not None:
        parent_terms = list(pathway_term.parents)
        parent_terms = [[parent_terms[i].name, parent_terms[i].id] for i in range(len(parent_terms))]
        for parent_id, parent_term in parent_terms:
            parent_path_lst.append([parent_id, parent_term, og_id, term_name])
            parent_path_lst = parse_parents(go, parent_id, og_id=og_id, parent_path_lst=parent_path_lst)
        if og_id is None:
            og_id = pathway_term.id
            og_id = ''
        outdf = pd.DataFrame(parent_path_lst, columns=['Parent_ID', 'Parent_name', 'OG_ID', 'OG_name'])
        return outdf
    else:
        print(f"No term with name '{pathway_name}' found in Gene Ontology")
        return pd.DataFrame()



pathways = 'your_pathway'
pathway_tree = pd.concat([parse_parents(go,p) for p in pathways]).drop_duplicates()

