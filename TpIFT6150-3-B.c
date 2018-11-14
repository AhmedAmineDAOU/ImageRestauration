/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  : DAOU Ahmed Amine & WARSHE Wrushabh         */
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

/*fonction utilisee dans le calcule de L'ISNR
  calcule les sommes des differences carrees des 
  matrices passées en parametres*/

float sum_diff(float** mat1, float** mat2, int length, int width) {

    float sum = 0.0;
    
    for( int i=0 ; i < length ; i++ ){
      for( int j=0 ; j < width ; j++ ) {
            sum += (mat1[i][j] - mat2[i][j])  * (mat1[i][j] - mat2[i][j]);
        }
    }

    return sum;
}

void clean(float** matRf,int length,int width){
  for (int i = 0 ; i < length ; i++){      
  for (int j = 0 ; j < width ; j++) {

    if (matRf[i][j]<0) matRf[i][j]=0.0;
    else if (matRf[i][j]>255) matRf[i][j]=255.0;

      }
    }
  }

void denoise(float** haar,int nbLevels,float threshold, int length, int width){

    for(int i = 0 ; i < length ; i++ )
        for(int  j=0 ; j < width ; j++ ) {
            if(i > length / powf(2, nbLevels) || j > length / powf(2, nbLevels)){
            
            if(fabs(haar[i][j]) < threshold)
                haar[i][j] = 0;
            else
                haar[i][j] = 255 * haar[i][j]/fabs(haar[i][j]) * (fabs(haar[i][j]) - threshold) / (255 - threshold);
        }}

}



int main(int argc,char** argv){
    int i,j,k,l;
    int length,width;
    int nbLevels;
    
    float threshold,ISNR,var;
       
    printf("Entrez la variance du bruit: ");
    scanf("%f",&var);
    
    printf("Entrez le nombre de niveaux a traiter : ");
    scanf("%d",&nbLevels);
        
    printf("Entrez le seuil : ");
    scanf("%f",&threshold);

    /* ouvrir l'image d'entree a garder pour l'isnr*/
    float** photograph = LoadImagePgm(NAME_IMG_IN, &length, &width);
    /*allocation memoires des autres matrices*/
    float** imageR = LoadImagePgm(NAME_IMG_IN, &length, &width);
    float** haar = fmatrix_allocate_2d(length, width);
    float** tmp1 = fmatrix_allocate_2d(length, width);
    
    /* ajouter du bruit a l'image d'entree (add_gaussian_noise) */
    add_gaussian_noise(imageR, length, width, var);
    

    /*transformee de Haar de imageR dans haar*/
    haar2D_complete(imageR, haar, nbLevels, length, width);


    /* debruiter l'image en seuillant les coefficients de Haar */
    denoise(haar, nbLevels, threshold,length, width);


    // Haar inverse de haar dans tmp1
    ihaar2D_complete(haar,tmp1,nbLevels, length,width);

  
    /* afficher l'ISBN */
    ISNR = 10 * log10(sum_diff(photograph, imageR, length, width) / sum_diff(photograph, tmp1, length, width));

    printf("fin-ISNR : %f\n", ISNR);
    
    /* sauvegarder les images */
    clean(tmp1,length,width);
    SaveImagePgm(NAME_IMG_OUT1, imageR, length, width);

    SaveImagePgm(NAME_IMG_OUT2, tmp1, length, width);

    /*liberation memoire*/
    free_fmatrix_2d(photograph);
    free_fmatrix_2d(imageR);
    free_fmatrix_2d(haar);
    free_fmatrix_2d(tmp1);
    
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");
    return 0;    
}






