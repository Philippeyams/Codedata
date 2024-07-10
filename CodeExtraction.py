import os
import numpy as np
import pandas as pd
import glob

chemin = "C:/Users/Hp/Desktop/Dossier Stage HSM/Post traitement DassFlow/Input brute/*.dat"
fic = glob.glob(chemin)

# Dossier de sortie
output_dir = "C:/Users/Hp/Desktop/Dossier Stage HSM/Post traitement DassFlow/Input"
os.makedirs(output_dir, exist_ok=True)  # Créer le dossier de sortie s'il n'existe pas

dataframes = {}

# Lire chaque fichier et le stocker dans le dictionnaire
for fichier in fic:
    df = pd.read_csv(fichier, delim_whitespace=True, skipinitialspace=True, header=1)
    dataframes[fichier] = df  
    print(f"Données du fichier {fichier}:")
    # print(df)

    columns_to_display = ["i", "x", "bathy"]
    if all(column in df.columns for column in columns_to_display):
        df_filtered = df[columns_to_display]
        print(df_filtered)

        base_filename = os.path.basename(fichier).replace('.dat', '_res.txt')
        output_file = os.path.join(output_dir, base_filename)
        
        # Enregistrer df_filtered en fichier texte sans les noms des colonnes
        df_filtered.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"DataFrame enregistré en {output_file} sans les noms des colonnes.")
    else:
        print(f"Les colonnes {columns_to_display} ne sont pas toutes présentes dans le fichier {fichier}.")

df1 = dataframes[fic[0]]
print("Données du premier fichier (accès via dictionnaire):")
#print(df1)
