/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  : DAOU Ahmed Amine & Warshe Wrudhabh         */
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

#define NAME_IMG_OUT1 "photograph_original_A"
#define NAME_IMG_OUT2 "photograph_degraded_withoutNoise_A" 
#define NAME_IMG_OUT3 "photograph_restored_withoutNoise_A"  
#define NAME_IMG_OUT4 "photograph_degraded_withNoise_A" 
#define NAME_IMG_OUT5 "photograph_restored_withNoise_A"  

/*applique un flou (filtre pass-bas (ilterR, filterI) de taille 
size_filter sur l'image (imageR,  imageI) et enregistre le resultat
dans l'image (gR,gI)*/
void flouer(float** gR , float** gI, float** imageR , float** imageI,float** filterR,float** filterI,int size_filter,int length,int width){

  for (int i = 0; i < length; i++) {

    for (int j = 0; j < width; j++) {

      if ( (i < size_filter/2.0 || i > length- size_filter/2.0) && (j < size_filter/2.0 || j > width- size_filter/2.0)){

        imageI[i][j] = filterI[i][j] = gR[i][j] = gI[i][j] = 0.0;
        filterR[i][j]= 1.0 / (size_filter*size_filter);
      }
      else imageI[i][j] = filterI[i][j] = gR[i][j] = gI[i][j] =filterR[i][j]= 0.0;
    }

  }
  /*multiplication de image * filtre dans le domaine frequenciel*/
  FFTDD(imageR, imageI, length, width);
  FFTDD(filterR, filterI, length, width);
  MultMatrix(gR, gI, imageR, imageI, filterR, filterI, length, width);

  /*retour dans le domaine spacial*/
  IFFTDD(gR,gI,length,width);
  IFFTDD(imageR, imageI, length, width);
  

}
/*la restauration suivant cette methode de descente de gradient
  genere des valeurs inferieures a 0 et superieures a 255
  qui affecte le resultat lors de l'enregistrement de l'image dans 
  le format pgm , cette fonction met les valeur >255 a 255
  et a 0 les valeurs <0    */
void clean(float** matRf,int length,int width){
  for (int i = 0 ; i < length ; i++){      
  for (int j = 0 ; j < width ; j++) {

    if (matRf[i][j]<0) matRf[i][j]=0.0;
    else if (matRf[i][j]>255) matRf[i][j]=255.0;

      }
    }
  }

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
/* restaure l'image  (gR,gI) a laquelle on a appliqué le filtre h (filterR, filterI)
   le resultat est stocké dans (matRf, matIf)
    imageR est passee en parametre pour calculer l'ISNR.
    */
void  restaurer(float** matRf, float**matIf, float** gR, float** gI,float** imageR,float** filterR, float**filterI, int pas, int nb_iterations,int length,int width ){
    float numerateur = 0.0, ISNR = 0.0;
    /* result convolution h*f */
  float** CI = fmatrix_allocate_2d(length,width);         
  float** CR = fmatrix_allocate_2d(length,width);

    /*result g- h*f*/
  float** parentheseI = fmatrix_allocate_2d(length,width);         
  float** parentheseR = fmatrix_allocate_2d(length,width);

    /*result h*(g_h*f)*/
  float** allI = fmatrix_allocate_2d(length,width);         
  float** allR = fmatrix_allocate_2d(length,width);


  int iteration=0;
    /*initialisaton matrice*/
  for (int i = 0; i < length; i++){ 
    for (int j = 0; j < width; j++){ 
      matRf[i][j]=gR[i][j];  //f0=g

      CI[i][j]=0.0;
      CR[i][j]=0.0;

      parentheseI[i][j]=0.0;
      parentheseR[i][j]=0.0;

      allR[i][j]=0.0;
      allI[i][j]=0.0;
        }
      }
  numerateur = sum_diff(imageR, gR, length, width);

  while(iteration < nb_iterations){

    FFTDD(matRf,matIf,length,width);

    MultMatrix(CR,CI,filterR,filterI,matRf,matIf,length,width);

    IFFTDD(CR,CI,length,width);

      

    substract(parentheseR,gR,CR,length,width);
    substract(parentheseI,gI,CI,length,width);


    FFTDD(parentheseR,parentheseI,length,width);



    MultMatrix(allR,allI,filterR,filterI,parentheseR,parentheseI,length,width);

    IFFTDD(allR,allI,length,width);
    IFFTDD(parentheseR,parentheseI,length,width);
    IFFTDD(matRf,matIf,length,width);

    add(matRf,allR,matRf,length,width);
    clean(matRf,length,width);
    ISNR = 10 * log10(numerateur / sum_diff(imageR, matRf, length, width));
    printf("%i%s%f\n",iteration," - ISNR ",ISNR);
    iteration++;  
  }
  free_fmatrix_2d(allR);
  free_fmatrix_2d(allI);

  free_fmatrix_2d(CR);
  free_fmatrix_2d(CI);

  free_fmatrix_2d(parentheseR);
  free_fmatrix_2d(parentheseI);

}



int main(int argc,char** argv)
 {
  int nb_iterations;
  int length,width;
  int pas=1;
  float var;
  int size_filter; /* taille du filtre servant a ajouter du flou a l'image d'entree */

   


  /* lire l'image d'entree */ 
  float** imageR = LoadImagePgm(NAME_IMG_IN, &length, &width);        /*photographe reelle*/
  float** imageI = fmatrix_allocate_2d(length,width);        /*photographe imaginaire*/

  float** filterR = fmatrix_allocate_2d(length,width);       /*filtre passe-bas reel*/
  float** filterI = fmatrix_allocate_2d(length,width);       /*filtre passe-bas imaginaire*/

  float** gR = fmatrix_allocate_2d(length,width);           /*photographe * filtre  (floué) reel*/
  float** gI = fmatrix_allocate_2d(length,width);           /*photographe * filtre  (floué) imaginaire*/

  float** matRf = fmatrix_allocate_2d(length,width);        /*Matrice Reelle resultat des iterations LANDWEBER*/
  float** matIf = fmatrix_allocate_2d(length,width);       /*Matrice Imaginaire resultat des iterations LANDWEBER*/          

  float** matNR = fmatrix_allocate_2d(length,width);        /*Matrice flou+bruit Reelle restauree*/
  float** matNI = fmatrix_allocate_2d(length,width);        /*Matrice flou+bruit Imaginaire restauree*/ 



	printf("Entrez la largeur du filtre passe bas : ");
  scanf("%d",&size_filter);

  /* Ajouter du flou a l'image d'entree */
  flouer(gR , gI, imageR , imageI, filterR, filterI, size_filter,length, width);
  /*enregistrer l'image flouee*/
  SaveImagePgm(NAME_IMG_OUT2,gR,length,width);

  printf("\nEntrez le nombre d'iterations pour LANDWEBER: ");
  scanf("%d",&nb_iterations);
  /*restaurer l'image flouer*/
  restaurer(matRf, matIf, gR,  gI, imageR, filterR, filterI, pas, nb_iterations, length, width );
  /*enregistrer l'image degradee apres restauration*/
  SaveImagePgm(NAME_IMG_OUT3,matRf,length,width);


  /* Ajouter du bruit a l'image floue : g = g + bruit = image + flou + bruit (add_gaussian_noise) */
 
  
	printf("Entrez la variance du bruit : ");
  scanf("%f",&var);
  /*ajouter du bruit gaussien*/
  add_gaussian_noise(gR, length, width, var);

  /*sauvegarder image contenant du flou + bruit*/
  SaveImagePgm(NAME_IMG_OUT4,gR,length,width);

  printf("\nEntrez le nombre d'itérations pour LANDWEBER: ");
  scanf("%d",&nb_iterations);

  /*restaurer l'image apres flou+bruit*/
  restaurer(matNR, matNI, gR,  gI, imageR, filterR, filterI, pas, nb_iterations, length, width );

  
  /*Sauvegarde de l'image flou+bruit restauree */
  SaveImagePgm(NAME_IMG_OUT5,matNR,length,width);

  

  /*Liberation memoire*/
  free_fmatrix_2d(imageR);
  free_fmatrix_2d(imageI);

  free_fmatrix_2d(filterR);
  free_fmatrix_2d(filterI);

  free_fmatrix_2d(gR);
  free_fmatrix_2d(gI);

  free_fmatrix_2d(matRf);
  free_fmatrix_2d(matIf);

  free_fmatrix_2d(matNR);
  free_fmatrix_2d(matNI);
  
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
