import requests
from Bio import Entrez
import os

Entrez.email = "tutu@yahoo.fr"


def id_genome(espece: str = "Homo sapiens") -> str:
    """
    Récupère l'ID d'un assemblage d'un génome pour une espèce donnée.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    params = {
    "db": "assembly",
    "term": f"{espece}[Organism]",
    "retmax": 20,
    "retmode": "json"
    }

    rez = requests.get(url, params=params)

    return rez.json()["esearchresult"]["idlist"]


def information_genome(aid: str) -> dict:
    """
    Récupère les informations sur un génome à partir de son accession number.
    """
    handle = Entrez.esummary(db="assembly", id=aid)
    summary = Entrez.read(handle)
    handle.close()

    return summary["DocumentSummarySet"]["DocumentSummary"][0]


def parser_ncbi_assembly(data: dict, colonnes: list[str]) -> str:
    """
    Parse les données ncbi pour en faire une ligne de tableau.
    """
    texte = ""
    
    for colonne in colonnes:
        if colonne in ["Sex", "Isolate"]:
            texte += f"{data["Biosource"][colonne]};"
        else:
            texte += f"{data[colonne]};"

    return texte


def telechargement(url: str, filename: str) -> None:
    """
    Télécharge et enregistre un fichier.
    """
    response = requests.get(url, stream=True)
    
    if response.status_code == 200:
        with open(filename, "w") as f:
            f.write(response.text)
    
    return


def generation_dataset(espece: str = "Homo sapiens") -> None:
    """
    Télécharge les informations relatives à un génome.
    """
    aids = id_genome(espece)
    donnees = []
    colonnes = ["ContigN50", "ScaffoldN50", "Coverage", "Sex", "Isolate", "BioSampleId", "AsmReleaseDate_GenBank", "FtpPath_GenBank"]
    texte  = ";".join(colonnes) + "\n"

    for aid in aids:
        donnees.append(information_genome(aid))

    for data in donnees:
        texte += parser_ncbi_assembly(data, colonnes)[:-1] + "\n"
    
    with open(f"{espece.lower().replace(" ", "_")}.txt", "w") as fout:
        fout.write(texte[:-1])

    return


def generation_dataset_cds(espece: str = "Homo sapiens") -> None:
    """
    Télécharge le fichier contenant le génome d'un espèce et le décompresse.
    """
    url = information_genome(id_genome(espece)[0])["FtpPath_GenBank"].replace("ftp://", "https://")
    genbank_url = url + f"/{url.split("/")[-1]}" +"_genomic.gbff.gz"
    genbank_filename = f"{espece.lower().replace(" ", "_")}.gbff.gz"
    telechargement(genbank_url, genbank_filename)
    os.system(f"gunzip -f {genbank_filename}")

    return