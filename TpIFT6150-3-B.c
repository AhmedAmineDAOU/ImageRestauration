/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  : DAOU Ahmed Amine & Wrushabh Warshe         */
/* Date    :                                            */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo3.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/  
/*------------------------------------------------*/
#define NAME_IMG_IN  "photograph"

#define NAME_IMG_OUT1 "photograph_bruite_B"
#define NAME_IMG_OUT2 "photograph_debruite_B"

void clean(float** matRf,int length,int width){
  for (int i = 0 ; i < length ; i++){      
  for (int j = 0 ; j < width ; j++) {

    if (matRf[i][j]<0) matRf[i][j]=0.0;
    else if (matRf[i][j]>255) matRf[i][j]=255.0;

      }
    }
  }

float sum_diff(float** mat1, float** mat2, int length, int width) {

    float sum = 0.0;
    
    for( int i=0 ; i < length ; i++ ){
      for( int j=0 ; j < width ; j++ ) {
            sum += (mat1[i][j] - mat2[i][j])  * (mat1[i][j] - mat2[i][j]);
        }
    }

    return sum;
}

int main(int argc,char** argv){
  	int i,j,k,l;
  	int length,width;
  	int nbLevels;
	
	float threshold, ISNR;
	float var;   
	float** haar;
	float** haar_inverse;
	float** imageR ;
	float** photographe;
	float** tmp1;

	
	printf("Entrez la variance du bruit: ");
  	scanf("%f",&var);
	
	printf("Entrez le nombre de niveaux a traiter N : ");
  	scanf("%d",&nbLevels);
		
	printf("Entrez le seuil TAU : ");
  	scanf("%f",&threshold);

	/* ouvrir l'image d'entree */
  	imageR = LoadImagePgm(NAME_IMG_IN, &length, &width);  /*photographe reelle a modifier*/
  	photographe = LoadImagePgm(NAME_IMG_IN, &length, &width);  /*photographe reelle a garder pour calcul ISNR*/
  	haar = fmatrix_allocate_2d(length,width);
	haar_inverse = fmatrix_allocate_2d(length,width);
	tmp1=fmatrix_allocate_2d(length,width);

	for(int i=0 ; i < length ; i++)  
    for(int j=0 ; j < width ; j++ ) 
 		haar_inverse[i][j]=0.0;
	

  	
	/* 1) etape ajouter du bruit a l'image d'entree (add_gaussian_noise) */

  	add_gaussian_noise(imageR, length, width, var);
  	SaveImagePgm(NAME_IMG_OUT1,imageR,length,width);


	
 	/*2) debruiter l'image en seuillant les coefficients de Haar */
	haar2D_complete(imageR, haar, nbLevels, length, width);

	for(int i=0;i< length / powf(2, nbLevels) ;i++)
        for(int j=0;j<width / powf(2, nbLevels);j++) {

            if(fabs( haar[i][j] ) < threshold)
                haar[i][j] = 0;
            else
                haar[i][j] = 255 * haar[i][j] / fabs(haar[i][j] ) * ( fabs(haar[i][j]) - threshold) / (255 - threshold);
	}

 	ihaar2D_complete(haar,haar_inverse,nbLevels,length,width);
	/* afficher l'ISBN */
	//clean(haar_inverse,length,width);

	/*ISNR = 10 * log10( sum_diff(photographe , imageR, length, width) / sum_diff(photographe, haar_inverse, length, width));
    printf("%s%f\n","fin - ISNR ",ISNR);*/
 	
 	/* sauvegarder les images */
 
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
