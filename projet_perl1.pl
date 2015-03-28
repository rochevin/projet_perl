#! /usr/bin/perl -w
#Auteur du script : Vincent ROCHER
#But du script : Traduire une séquence ADN en protéine selon les 6 cadres de lecture
#Fonctionne avec du multifasta -> créer un fichier .fasta par séquence avec les 6 cadres

use strict;
use warnings;

#########################INPUT###############################

print "####################################################\n";
print "#                                                  #\n";
print "#   Program : Translate DNA sequence from fasta    #\n";
print "#    file into protein sequence with six frame     #\n";
print "#             Author : Vincent ROCHER              #\n";
print "#                                                  #\n";
print "####################################################\n";
print "Fasta or multi fasta file directory : ";
my $file_write=<>;
#Appel de la fonction fasta_file_write, enregistre la sortie dans une table de hachage
my (%SEQ)=fasta_file_write($file_write);
######TRADUCTION EN PROTEINE#########
#%SEQ contient toutes les séquences issues du fasta, avec leur nom comme clé et la séquence comme valeur
#On fait tourner un foreach tant qu'il y a des séquences dans la table
foreach my$elmt(keys %SEQ) {
	#Protein_list contiendra les 6 cadres de lecture de la séquence et son nom
	my@proteine_list = (); 
	push @proteine_list, $elmt; #On enregistre le nom de la séquence comme premier élément
	push @proteine_list, trad_peptide($SEQ{$elmt}, 0); #Traduit la protéine avec le premier cadre de lecture
	push @proteine_list, trad_peptide($SEQ{$elmt}, 1); #Traduit la protéine avec le second cadre de lecture
	push @proteine_list, trad_peptide($SEQ{$elmt}, 2); #Traduit la protéine avec le troisième cadre de lecture

	my$reverse = revcom($SEQ{$elmt}); #On appelle la fonction revcom qui va renvoyer le brin complémentaire de la séquence
	
	push @proteine_list, trad_peptide($reverse, 0); #Traduit la protéine avec le premier cadre de lecture du brin complémentaire
	push @proteine_list, trad_peptide($reverse, 1); #Traduit la protéine avec le second cadre de lecture du brin complémentaire
	push @proteine_list, trad_peptide($reverse, 2); #Traduit la protéine avec le troisième cadre de lecture du brin complémentaire
	

	#On donne des informations à l'utilisateur
	print "####################################################\n";
	fasta_print(@proteine_list); #On envoit la liste avec les 6 cadres à la fonction qui va permettre l'écriture dans un multi-fasta

}

######FONCTIONS########################
#But de la fonction : Ouvrir un fichier (multi)fasta
sub fasta_file_write {
	my($file_write_name)=@_; #On recupère le nom du fichier
	#On ouvre le fichier et on affiche un message d'erreur en cas d'echec
	unless ( open(file_write_name, $file_write_name) ) {
        print STDERR "\"$file_write_name\" can't be found...\n\n";
        exit;
    }
    #On enregistre le contenu dans une liste et on ferme
	my@file_write_seq= <file_write_name>;
	close file_write_name;
	#On déclare les variables vides au préalables
	my$data_name="";
	my$data_seq="";
	my%list_sequence=(); #Contiendra l'ensemble des séquences avec leurs nom, sera renvoyer en return
	#On parcours chaque ligne du fichier
	foreach my$line (@file_write_seq) {
		if ( $line =~ /^\s*$/ ){
            next; #On passe à la ligne suivante si la ligne ne contient que des espaces
        }
        #On vérifie si la ligne commence par '>'
        elsif ( $line =~ /^>/ ) {
        	#Si oui et que la variable contient la séquence, ça veut dire qu'on à affaire à une nouvelle séquence,
        	#on enregistre la précédente dans la table si elle existe
        	if ($data_seq) {
        		#On applique certaines expression pour nettoyer le contenu d'élements indésirables
                $data_seq =~ s/\s//g; #espaces quelquonques
                $data_seq =~ s/[0-9]//g; #chiffres
                $data_seq =~ tr/Uu/Tt/; #Si on a affaire à de l'ARN 
                #On enregistre dans la table avec comme clé le nom
                $list_sequence{$data_name}=$data_seq;
            }
            #On enregistre le nom de la séquence, qu'on utilisera comme clé
            $data_name = $line;
            $data_name =~ s/^.//; #On suprimme le premier caractère ">"
            chomp $data_name; #On enlève le \n en fin de ligne
            #On nettoie la variable pour ne pas qu'elle contienne l'ancienne séquence + la nouvelle
            $data_seq="";
            
        }
        #Sinon, on considère qu'on est sur une ligne avec de la séquence ADN
        else {
        	#On ajoute la ligne à la variable qui contiendra la séquence en entier
            $data_seq .= $line;
        }
	}
	#Ces 4 dernières lignes servent à enregistrer la dernière séquence :
    $data_seq =~ s/\s//g; 
    $data_seq =~ s/[0-9]//g;
    $data_seq =~ tr/Uu/Tt/;
    $list_sequence{$data_name}=$data_seq;
    #On renvoi la table avec l'ensemble des séquences et leurs noms
    return (%list_sequence);
}
#But de la fonction : Créer le brin complémentaire
sub revcom {
    my($adn)=@_; #On récupère
    my $revcom= reverse $adn; #On inverse le brin
    $revcom =~ tr/ACGTacgt/TGCAtgca/; #On change toutes les lettres
    #On renvoit la nouvelle variable qui contient la fonction
    return $revcom;
}
#But de la fonction : Recupérer la séquence avec son cadre de lecture, et envoyer codon par codon
#à la fonction qui contient le code génétique
sub trad_peptide {
    my($sequence, $debut)=@_; #On recupère la séquence et son cadre 
    #On determine la taille de la séquence pour construire la sous chaine
    my$fin= length($sequence);
    #Sous chaine qui contient la séquence commençant au bon cadre de lecture
    my$cadre=substr($sequence, $debut , $fin - $debut);
    
    #On déclare une variable vide qui va contenir la séquence protéique
    my$proteine="";
    #On fait une boucle qui va travailler sur la séquence et la diviser en triplet
    for (my$i=0; $i<(length($cadre) - 2); $i+=3) {
    		#On enregistre progressivement chaque nouvel acide aminé
            $proteine .= codon_into_aa(substr($cadre,$i,3));
        }
    #On renvoit la séquence protéique
    return $proteine;
}
#But de la fonction : Dictionnaire avec les triplets comme clé et l'acide aminé correspondant comme valeur
sub codon_into_aa {
    my($codon)=@_; #On récupère
	my%genetic_code = (
		# Tryptophan
	    'TGG' => 'W',  
	    # Methionine
	    'ATG' => 'M',  
	    # Aspartic Acid
	    'GAC' => 'D','GAT' => 'D', 
	    # Glutamic Acid   
	    'GAA' => 'E','GAG' => 'E', 
	    # Asparagine
	    'AAC' => 'N','AAT' => 'N',    
	    # Lysine
	    'AAA' => 'K','AAG' => 'K',
	    # Phenylalanine
	    'TTC' => 'F','TTT' => 'F',
	    # Tyrosine
	    'TAC' => 'Y','TAT' => 'Y',  
	    # Histidine
	    'CAC' => 'H','CAT' => 'H',   
	    # Glutamine 
	    'CAA' => 'Q','CAG' => 'Q',     
	    # Cysteine
	    'TGC' => 'C','TGT' => 'C',    
	    # Isoleucine
	    'ATA' => 'I','ATC' => 'I','ATT' => 'I',
	    # Stop    
	    'TAA' => '*','TAG' => '*','TGA' => '*',   
	    # Proline    
	    'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P',
	    # Threonine
	    'ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
	    # Valine
	    'GTA' => 'V','GTC' => 'V','GTG' => 'V','GTT' => 'V',
	    # Alanine
	    'GCA' => 'A','GCC' => 'A','GCG' => 'A','GCT' => 'A',
	    # Glycine
	    'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G',
	    # Serine
	    'TCA' => 'S','TCC' => 'S','TCG' => 'S','TCT' => 'S','AGC' => 'S','AGT' => 'S',
	    # Leucine
	    'CTA' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L','TTA' => 'L','TTG' => 'L',
	    # Arginine  
	    'CGA' => 'R','CGC' => 'R','CGG' => 'R','CGT' => 'R','AGA' => 'R','AGG' => 'R',   
	 );
	#On renvoit la valeur correspondant à la clé(triplet) envoyée par la fonction trad_peptide
	return $genetic_code{$codon};
}
#But de la fonction : Écrire dans un fichier multifasta les protéines correspondants au 6 cadres de lecture de chaque séquence
sub fasta_print {  
        
    my(@list_prot)=@_;
    #On recupère le nom de la séquence et on l'enlève de la liste
    my$name= shift @list_prot;
    $name =~ s/\s//g;
    #On créer un fichier en mode écriture
    unless ( open(file_write, ">$name.fasta" ) ) {
        print STDERR "$name.fasta can't be created\n" ;
        exit;
    }
    #On utilise une double boucle for
    #La première va faire un tour autant de fois qu'il y a de cadres, soit 6
    for (my $i = 0; $i < 6; $i++) {
        print file_write ">".$name."|Frame=".($i+1)."\n"; #Et écris le nom de la séquence ainsi que son cadre
        #La deuxième boucle va diviser la séquence protéique en morceau de 70 Acides aminés
        for ( my $pos = 0 ; $pos < length($list_prot[$i]) ; $pos +=70) {
        	my$print = substr($list_prot[$i], $pos, 70); 
        	print file_write "$print\n";
        }
    }
    #Envoi le nom du fichier produit à l'utilisateur via le terminal
    print "#File : $name.fasta\n";
    #Ferme le fichier une fois l'opération terminée
	close file_write;

}
#Fin du programme