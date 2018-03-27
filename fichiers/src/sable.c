#include "global.h"
#include "compute.h"
#include "graphics.h"
#include "debug.h"
#include "ocl.h"
#include "scheduler.h"

#include <stdbool.h>

static long unsigned int *TABLE=NULL;

static volatile int changement;
static volatile int changement_tile;

static unsigned long int max_grains;

#define table(i,j) TABLE[(i)*DIM+(j)]

#define RGB(r,v,b) (((r)<<24|(v)<<16|(b)<<8))

void sable_init()
{
  TABLE=calloc(DIM*DIM,sizeof(long unsigned int));
}

void sable_finalize()
 {
   free(TABLE);
 }


///////////////////////////// Production d'une image
void sable_refresh_img()
{
  unsigned long int max = 0;
    for (int i = 1; i < DIM-1; i++)
      for (int j = 1; j < DIM-1; j++)
	{
	  int g =table(i,j);
	  int r,v,b;
	  r = v = b = 0;
	    if ( g == 1)
	      v=255;
	    else if (g == 2)
	      b = 255;
	    else if (g == 3)
	      r = 255;
	    else if (g == 4)
	      r = v = b = 255 ;
	    else if (g > 4)
	      r = b = 255 - (240 * ((double) g) / (double) max_grains);

	    cur_img (i,j) = RGB(r,v,b);
	    if (g > max)
	      max = g;
	}
    max_grains = max;
}



///////////////////////////// Configurations initiales

void sable_draw_4partout(void);

void sable_draw (char *param)
{
  char func_name [1024];
  void (*f)(void) = NULL;

  sprintf (func_name, "draw_%s", param);
  f = dlsym (DLSYM_FLAG, func_name);

  if (f == NULL) {
    printf ("Cannot resolve draw function: %s\n", func_name);
    f = sable_draw_4partout;
  }

  f ();
}


void sable_draw_4partout(void){
  max_grains = 8;
  for (int i=1; i < DIM-1; i++)
    for(int j=1; j < DIM-1; j++)
      table (i, j) = 4;
}


void sable_draw_DIM(void){
  max_grains = DIM;
   for (int i=DIM/4; i < DIM-1; i+=DIM/4)
     for(int j=DIM/4; j < DIM-1; j+=DIM/4)
       table (i, j) = i*j/4;
}


 void sable_draw_alea(void){
   max_grains = DIM;
 for (int i= 0; i < DIM>>2; i++)
   {
     table (1+random() % (DIM-2) , 1+ random() % (DIM-2)) = random() % (DIM);
   }
}


///////////////////////////// Version séquentielle simple (seq)

static inline void compute_new_state (int y, int x)
{
  if (table(y,x) >= 4)
    {
      unsigned long int div4 = table(y,x) / 4;
      table(y,x-1)+=div4;
      table(y,x+1)+=div4;
      table(y-1,x)+=div4;
      table(y+1,x)+=div4;
      table(y,x)%=4;
      changement = 1;
      changement_tile = 1;
    }
}

static void traiter_tuile (int i_d, int j_d, int i_f, int j_f)
{
  PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

  for (int i = i_d; i <= i_f; i++)
    for (int j = j_d; j <= j_f; j++)
      compute_new_state (i, j);
}

// Renvoie le nombre d'itérations effectuées avant stabilisation, ou 0
unsigned sable_compute_seq (unsigned nb_iter)
{

  for (unsigned it = 1; it <= nb_iter; it ++) {
    changement = 0;
    // On traite toute l'image en un coup (oui, c'est une grosse tuile)
    traiter_tuile (1, 1, DIM - 2, DIM - 2);
    if(changement == 0)
      return it;
  }
  return 0;
}

/* Déroulage de boucle en 2 */
static void traiter_tuile2 (int i_d, int j_d, int i_f, int j_f) {

    PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

    int modj = (j_f - (j_d + 1)) % 2;

    for (int i = i_d; i <= i_f; i++) {
        for (int j = j_d; j <= j_f - modj; j+=2) {
            compute_new_state (i, j);
            compute_new_state (i, j+1);
        }
        if (modj ==1 ) compute_new_state (i, j_f);
    }
}

unsigned sable_compute_seq2 (unsigned nb_iter)
{

  for (unsigned it = 1; it <= nb_iter; it+=2) {
    changement = 0;
    traiter_tuile2 (1, 1, DIM - 2, DIM - 2);
    if(changement == 0)
      return it;
  }
  return 0;
}

/* Déroulage de boucle en 3 */
static void traiter_tuile3 (int i_d, int j_d, int i_f, int j_f) {

    PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

    int mod = (j_f - j_d + 1) % 3;
    for (int i = i_d; i <= i_f; i++) {
        for (int j = j_d; j <= j_f - mod; j+=3) {
            compute_new_state (i, j);
            compute_new_state (i, j+1);
            compute_new_state (i, j+2);
        }
        if (mod != 0) {
            compute_new_state (i, j_f);
            compute_new_state (i, j_f - 1);
        }
    }
}

unsigned sable_compute_seq3 (unsigned nb_iter)
{

  for (unsigned it = 1; it <= nb_iter; it+=3) {
    changement = 0;
    traiter_tuile3 (1, 1, DIM - 2, DIM - 2);
    if(changement == 0)
      return it;
  }
  return 0;
}

/* Déroulage de boucle en 4 */
static void traiter_tuile4 (int i_d, int j_d, int i_f, int j_f) {

    PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

    int mod = (j_f - j_d + 1) % 4;
    for (int i = i_d; i <= i_f; i++) {
        for (int j = j_d; j <= j_f - mod; j+=4) {
            compute_new_state (i, j);
            compute_new_state (i, j+1);
            compute_new_state (i, j+2);
            compute_new_state (i, j+3);
        }
        if (mod != 0) {
            compute_new_state (i, j_f);
            compute_new_state (i, j_f - 1);
            compute_new_state (i, j_f - 2);
        }
    }
}

unsigned sable_compute_seq4 (unsigned nb_iter)
{

  for (unsigned it = 1; it <= nb_iter; it+=3) {
    changement = 0;
    traiter_tuile4 (1, 1, DIM - 2, DIM - 2);
    if(changement == 0)
      return it;
  }
  return 0;
}

/* Déroulage de boucle via GCC*/

#pragma GCC push_options
#pragma GCC optimize("unroll-all-loops")
static void traiter_tuilegcc (int i_d, int j_d, int i_f, int j_f) {

    PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

    for (int i = i_d; i <= i_f; i++) {
        for (int j = j_d; j <= j_f; j++)
            compute_new_state (i, j);
    }
}
#pragma GCC pop_options

unsigned sable_compute_seqgcc (unsigned nb_iter)
{

  for (unsigned it = 1; it <= nb_iter; it+=3) {
    changement = 0;
    traiter_tuilegcc (1, 1, DIM - 2, DIM - 2);
    if(changement == 0)
      return it;
  }
  return 0;
}




///////////////////////////// Version séquentielle tuilée (tiled)


static unsigned tranche = 0;

unsigned sable_compute_tiled (unsigned nb_iter)
{
  tranche = DIM / GRAIN;

  for (unsigned it = 1; it <= nb_iter; it ++) {

    // On itére sur les coordonnées des tuiles
    for (int i=0; i < GRAIN; i++)
      for (int j=0; j < GRAIN; j++)
	{
	  traiter_tuile (i == 0 ? 1 : (i * tranche) /* i debut */,
			 j == 0 ? 1 : (j * tranche) /* j debut */,
			 (i + 1) * tranche - 1 - (i == GRAIN-1)/* i fin */,
			 (j + 1) * tranche - 1 - (j == GRAIN-1)/* j fin */);
	}
  }

  return 0;
}

/* Version récursive tuilée */

unsigned int MAX_DEPTH = 32;

static void traiter_tuile_rec(int i, int j, int tranche, int depth) {
    if (i<0 || i>=GRAIN || j<0 || j>=GRAIN) return;

    int tmp = 0;
    int i_d = i == 0 ? 1 : (i * tranche); /* i debut */
    int i_f = (i + 1) * tranche - 1 - (i == GRAIN-1); /* i fin */
    int j_d = j == 0 ? 1 : (j * tranche); /* j debut */
    int j_f = (j + 1) * tranche - 1 - (j == GRAIN-1);/* j fin */
    PRINT_DEBUG ('c', "tuile [%d-%d][%d-%d] traitée\n", i_d, i_f, j_d, j_f);

    do {
        changement = 0;
        traiter_tuile2(i_d, j_d, i_f, j_f);
        if (changement == 1) tmp = 1;
    } while (changement == 1);

    if(tmp == 1 && depth < MAX_DEPTH){
    traiter_tuile_rec(i, j+1, tranche, ++depth);
    traiter_tuile_rec(i, j-1, tranche, ++depth);
    traiter_tuile_rec(i-1, j, tranche, ++depth);
    traiter_tuile_rec(i+1, j, tranche, ++depth);
  }
}

unsigned sable_compute_rec_tiled(unsigned nb_iter) {
    tranche = DIM / GRAIN;
    int depth = 0;

    for (unsigned it = 1; it <= nb_iter; it ++) {
        changement_tile = 0;
        for (int i=0; i < GRAIN; i++) {
            for (int j=0; j < GRAIN; j++) {
                traiter_tuile_rec(i, j, tranche, depth);
            }
        }
        if(changement_tile == 0) return it;
    }
    return 0;
}

/* Version parallèle lignes omp task */
unsigned sable_compute_omptask_line(unsigned nb_iter) {
    unsigned tmp = 0;
    tranche = DIM / GRAIN;
    int i = 0;
    #pragma omp parallel
    #pragma omp single
    for (unsigned it = 1; it <= nb_iter && tmp == 0; it++) {
        changement = 0;
        for (i = 0; i < GRAIN; i += 3) {
            #pragma omp task firstprivate(i, changement) depend(out : table(i, 0))
            {
                for (int j = 0; j < GRAIN; j++) {
                    traiter_tuile2(i == 0 ? 1 : (i * tranche) /* i debut */,
                                   j == 0 ? 1 : (j * tranche) /* j debut */,
                                  (i + 1) * tranche - 1 - (i == GRAIN - 1) /* i fin */,
                                  (j + 1) * tranche - 1 - (j == GRAIN - 1) /* j fin */);
                }
            }
        }
    for (int k = 1; k < 3; k++) {
        for (i = k; i < GRAIN; i += 3) {
            #pragma omp task firstprivate(i, changement) depend(in : table(i - 1, 0), table(i + 2, 0)) depend(out : table(i, 0))
            {
                for (int j = 0; j < GRAIN; j++) {
                    traiter_tuile2(i == 0 ? 1 : (i * tranche) /* i debut */,
                                   j == 0 ? 1 : (j * tranche) /* j debut */,
                                  (i + 1) * tranche - 1 - (i == GRAIN - 1) /* i fin */,
                                  (j + 1) * tranche - 1 - (j == GRAIN - 1) /* j fin */);
                }
            }
        }
    }
    #pragma omp taskwait
    if (changement == 0) tmp = it;
    }
    return tmp;
}

/* Version parallele en ligne omp for
 * Découpage en ligne :
 * 111
 * 222
 * 333
 */
unsigned sable_compute_ompfor_line (unsigned nb_iter) {

    for (unsigned it = 1; it <= nb_iter; it++) {
        changement = 0;

        #pragma omp parallel
        #pragma omp for
            for (int j = 1; j <= DIM - 2; j+=3) {
                traiter_tuile2(j, 1, j, DIM - 2);
            }

        #pragma omp for
            for (int j = 2; j <= DIM - 2; j+=3) {
                traiter_tuile2(j, 1, j, DIM - 2);
            }

        #pragma omp for
            for (int j = 3; j <= DIM - 2; j+=3) {
                traiter_tuile2(j, 1, j, DIM - 2);
            }

        if (changement == 0) return it;
        }
    return 0;
}

/* Version parallele tuilee omp for
 * Découpage en :
 * 123
 * 456
 * 789
 */
unsigned sable_compute_ompfor_tiled33 (unsigned nb_iter) {
    tranche = DIM / GRAIN;
    unsigned tmp = 0; //Contourne GCC
    #pragma omp parallel
    for (unsigned it = 1; it <= nb_iter && tmp == 0; it++) {
        changement = 0;
        for (int modi = 0; modi < 3; modi++) {
            for (int modj = 0; modj < 3; modj++) {
            #pragma omp for collapse(2)
                for (int i = modi; i < GRAIN; i+=3) {
                    for (int j = modj; j < GRAIN; j+=3) {

                        traiter_tuile2 (i == 0 ? 1 : (i * tranche) /* i debut */,
                                        j == 0 ? 1 : (j * tranche) /* j debut */,
                                       (i + 1) * tranche - 1 - (i == GRAIN-1)/* i fin */,
                                       (j + 1) * tranche - 1 - (j == GRAIN-1)/* j fin */);
                    }
                }
            }
        }
        if (changement == 0) tmp = it;
    }
    return tmp;
}

/*
 * 1234
 * 5678
 * 3412
 * 7856
 */

unsigned sable_compute_ompfor_tiled44(unsigned nb_iter)
{
  tranche = DIM / GRAIN;

  for (unsigned it = 1; it <= nb_iter; it ++) {
    changement = 0;
    #pragma omp parallel
    for (int modi = 0; modi < 2; modi++){
      for (int modj = 0; modj < 4; modj++){
        #pragma omp for schedule(static)
        // On itére sur les coordonnées des tuiles
        for (int i=modi; i < GRAIN; i+=2)
        for (int j=(modj+i%4)%4; j < GRAIN; j+=4){
          traiter_tuile (i == 0 ? 1 : (i * tranche) /* i debut */,
          j == 0 ? 1 : (j * tranche) /* j debut */,
          (i + 1) * tranche - 1 - (i == GRAIN-1)/* i fin */,
          (j + 1) * tranche - 1 - (j == GRAIN-1)/* j fin */);
        }
      }
    }
    if (changement == 0){
      return it;
    }
  }

  return 0;
}

///////////////////////////// Version utilisant un ordonnanceur maison (sched)

unsigned P;

void sable_init_sched ()
{
  sable_init();
  P = scheduler_init (-1);
}

void sable_finalize_sched ()
{
  sable_finalize();
  scheduler_finalize ();
}

static inline void *pack (int i, int j)
{
  uint64_t x = (uint64_t)i << 32 | j;
  return (void *)x;
}

static inline void unpack (void *a, int *i, int *j)
{
  *i = (uint64_t)a >> 32;
  *j = (uint64_t)a & 0xFFFFFFFF;
}

static inline unsigned cpu (int i, int j)
{
  return 1;
}

static inline void create_task (task_func_t t, int i, int j)
{
  scheduler_create_task (t, pack (i, j), cpu (i, j));
}

//////// First Touch

static void zero_seq (int i_d, int j_d, int i_f, int j_f)
{

  for (int i = i_d; i <= i_f; i++)
    for (int j = j_d; j <= j_f; j++)
      next_img (i, j) = cur_img (i, j) = 0 ;
}

static void first_touch_task (void *p, unsigned proc)
{
  int i, j;

  unpack (p, &i, &j);

  //PRINT_DEBUG ('s', "First-touch Task is running on tile (%d, %d) over cpu #%d\n", i, j, proc);
  zero_seq (i * tranche, j * tranche, (i + 1) * tranche - 1, (j + 1) * tranche - 1);
}

void sable_ft_sched (void)
{
  tranche = DIM / GRAIN;

  for (int i = 0; i < GRAIN; i++)
    for (int j = 0; j < GRAIN; j++)
      create_task (first_touch_task, i, j);

  scheduler_task_wait ();
}

//////// Compute

static void compute_task (void *p, unsigned proc)
{
  int i, j;

  unpack (p, &i, &j);

  //PRINT_DEBUG ('s', "Compute Task is running on tile (%d, %d) over cpu #%d\n", i, j, proc);
	  traiter_tuile (i == 0 ? 1 : (i * tranche) /* i debut */,
			 j == 0 ? 1 : (j * tranche) /* j debut */,
			 (i + 1) * tranche - 1 - (i == GRAIN-1)/* i fin */,
			 (j + 1) * tranche - 1 - (j == GRAIN-1)/* j fin */);
}

unsigned sable_compute_sched (unsigned nb_iter)
{
  tranche = DIM / GRAIN;

  for (unsigned it = 1; it <= nb_iter; it ++) {
    changement=0;
    for (int i = 0; i < GRAIN; i++)
      for (int j = 0; j < GRAIN; j++)
	create_task (compute_task, i, j);

    scheduler_task_wait ();

    if (changement == 0)
      return it;
  }

  return 0;
}
