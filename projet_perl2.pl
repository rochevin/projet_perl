#! /usr/bin/perl -w
#Auteur du script : Vincent ROCHER
#But du script : Extraire des informations d'un fichier .vcf pour faire un graphique sous R


use strict;
use warnings;
use Statistics::R;
use File::Temp qw/ tempfile /;
######CORPS DU CODE######
print "####################################################\n";
print "#                                                  #\n";
print "#     Program : Parsing variant calling format     #\n";
print "#     file to extract SNV and generate pie for     #\n";
print "#            Total and all chromosomes             #\n";
print "#             Author : Vincent ROCHER              #\n";
print "#                                                  #\n";
print "####################################################\n";
print "Enter your .vcf format file directory : ";
my$fichier=<>; #On demande à l'utilisateur de nommer le fichier
chomp $fichier; #On chomp le \n 
my%each_chms=open_vcf($fichier); #On envoi le nom du fichier à la fonction open_vcf et on récupère un double tableau associatif

#####FICHIER R TEMPORAIRE#####

my ($fh, $filename) = tempfile(SUFFIX => '.R'); #On créer un fichier temporaire avec l'extension .R
my @nom_fichier = split(/\./, $fichier); #Suprimme le .vcf du nom du fichier pour l'affichage en .pdf -> Evite le example.vcf.pdf
#On commence a écrire dans le fichier temporaire :
print $fh "pdf(\"$nom_fichier[0].pdf\")\n"; 
foreach my$chms(sort {substr($a, 3, 2) <=> substr($b, 3, 2)} keys %each_chms) { #Parcours du tableau associatif
	my@value_forR=(); #Liste vide qui va contenir les SNV de chaque chromosome à chaque tour du foreach
	my@value_name_forR=(); #Liste vide qui va contenir le nom des SNV de chaque chromosome à chaque tour du foreach
	my@colorsR=(); #Liste vide qui va contenir les couleurs correspondant aux SNV
	my $name = ""; #Variable vide qui contiendra le titre du graphique
   	if ($chms eq "chr23") { #On regarde si la clé a la valeur chr23,chr24 ou chr25 pour donner les bons titres dans le pdf
        			$name= "chrX";
    }
    elsif ($chms eq "chr24") {
        			$name= "chrY";
    }
    elsif ($chms eq "chr25") {
        			$name= "Total";
    }
    else {
    	$name = $chms; #Sinon c'est que la clé a le bon nom et qu'on peut l'afficher tel quel dans le pdf
    }
	print ("\n$name-> | "); #La on print dans le terminal le chromosome...
	#Puis un deuxième foreach pour parcourir les éléments du deuxième tableau associatif contenu dans le premier
	foreach (sort keys %{$each_chms{$chms}}) {
		print ("$_ : $each_chms{$chms}{$_} | "); #On print dans le terminal chaque clé (SNV) et leur valeur correspondante
		#La on va inscrire dans les listes les éléments du double tableau associatif pour l'envoyer au script R
		push @value_forR, $each_chms{$chms}{$_}; #On push chaque élément du tableau dans la liste correspondant aux valeurs
		push @value_name_forR, "\"".$_."\""; #On push chaque clé du tableau (SNV) dans la liste correspondant aux noms
		push @colorsR, "\"".SNV_into_color($_)."\""; #On push la couleur correspondant au SNV grâce à la fonction SNV_into_colors
		
    }
    print ("\n");
    #La on prépare les vecteurs pour R -> on join la liste dans le vecteur séparé par une virgule
	my $loadvect1="vect <- c (". join(",", @value_forR) .")" ; #Le vecteur contient la liste contenant le nombre de SNV dans chaque chromosome
	my $loadvect2="lbls <- c (". join(",", @value_name_forR) .")" ; #Le vecteur contient la liste contenant les noms des SNV dans chaque chromosome
	my $loadvect3="lblscolors <- c (". join(",", @colorsR) .")" ; #Le vecteur contient les couleurs correspondants aux différents SNV -> sert à garder les mêmes couleurs peu-importe le graphique
	#La on print les 3 vecteurs dans le fichier
	print $fh "$loadvect1\n"; #vect
	print $fh "$loadvect2\n"; #lbls
	print $fh "$loadvect3\n"; #lblscolors
	#On effectue un simple calcul pour afficher le pourcentage de chaque élément
	print $fh "pct <- round(vect/sum(vect)*100)\n";
	#Et on le colle aux valeurs pour l'afficher à coté
	print $fh "lbls <- paste(lbls, pct)\n";
	print $fh "lbls <- paste(lbls,\"%\",sep=\"\")\n"; #On met un % pour l'affichage 
	print $fh "pie(vect, labels = lbls, main=\"$name\", col = lblscolors)\n"; #la on fait le graphique en utilisant les 3 vecteurs
}
print $fh "dev.off()\n"; #On fini d'écrire dans le fichier et on indique qu'on a terminé de faire des graphiques avec un dev.off() dans R

########INTERACTION AVEC R#########
#on déclare une nouvelle session R
my $R = Statistics::R->new();
$R->startR;
$R->run_from_file($filename); #On lis le fichier temporaire avec R et on l'éxecute
$R->stopR(); #On stop la session R

##########FONCTIONS###########

sub open_vcf {
	my($file_name)=@_; #On recupère le nom du fichier

	unless ( open(VCF_CONTENT, $file_name) ) {
        print STDERR "Impossible d'ouvrir le fichier \"$file_name\"\n\n";
        exit;
    }

	my@file_vcf= <VCF_CONTENT>; #On enregistre le fichier dans une liste
	close VCF_CONTENT; #On ferme le fichier

	my%chms_content=(); #On déclare un tableau associatif vide
	foreach my$line (@file_vcf) { #Pour chaque ligne du fichier
		if ( $line =~ /^\s*$/ ){
            next; #Si la ligne est vide, on next
        }
        elsif ( $line =~ /^#/ ) {
            next; #Sinon si elle contient # comme premier caractère, on next
        }
        else { #Sinon
        	my @line_content = split("\t", $line); #On eclate chaque élément de la ligne séparé par une tabulation dans une liste
        	if ((length($line_content[0])==1)||(length($line_content[0])==2)) {
        		$line_content[0] = "chr".$line_content[0]; #Au cas ou "chr" n'est pas présent et qu'il n'y a que les numeros, on l'ajoute
        	}
        	if (($line_content[3] =~ /^[A-Za-z]$/)&&($line_content[4] =~ /^[A-Za-z]$/)) { #Si l'élement 3 et 4 de la ligne sont une simple lettre
        		my$local_SNV = uc($line_content[3].$line_content[4]); #On enregistre l'élement 3 et 4 ensemble dans une variable
        		if ($line_content[0] eq "chrX") { #Si l'élement 0 de la ligne est chrX, on le remplace par chr23 -> va servir pour le tri de la liste et permettre d'afficher
        			$line_content[0]= "chr23"; #les chromosomes dans l'ordre dans le pdf
        		}
        		if ($line_content[0] eq "chrY") { #Pareil pour chrY
        			$line_content[0]= "chr24";
        		}
        		$chms_content{$line_content[0]}{$local_SNV}++; #C'est là ou on compte le nombre de SNV dans chaque chromosome, qu'on enregistre dans une double table de hachage
        		$chms_content{"chr25"}{$local_SNV}++; #La on enregistre à chaque fois dans chr25 qui est le total, pour compter le nombre total de SNV
        	} 
        }
        
	}
	
    return (%chms_content); #On renvoit le double tableau associatif
}
#Dictionnaire de couleurs qui va servir à ce que chaque SNV ai la même couleur peu importe le graphique
sub SNV_into_color {
    my($SNV)=@_; #On récupère
	my%color = (
	    'AT' => 'red',
	    'AG' => 'red3',
	    'AC' => 'red4',
	    'TG' => 'purple', 
	    'TC' => 'purple3', 
	    'TA' => 'purple4', 
	    'GT' => 'blue', 
	    'GC' => 'blue3', 
	    'GA' => 'blue4', 
	    'CT' => 'green', 
	    'CG' => 'green3', 
	    'CA' => 'green4',  
	 );
	#On renvoit la couleur correspondant au SNV 
	return $color{$SNV}; 
}
#Fin du programme