/* Программе в аргументах командной строки задается Р -кол-во процессов. 
 Главный процесс должен создать Р процессов, которые работают параллельно 
с разделяемой памятью как переменной  и печатают свой номер последовательно  
и результат вычислений согласно нумерации процессов. 
Родитель дожидается завершения всех созданных им процессов, выводит на 
sddout число 0 и сам завершает свою работу с кодом завершения 0. Для 
синхронизации процессов использовать семафоры. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/shm.h>          //Добавлена shared memory
#include <sys/sem.h>


/*
struct sembuf {
short sem_num;   // номер семафора в векторе
short sem_op;    // производимая операция
short sem_flg;   // флаги операции

union semun {
int val; 		// значение одного семафора
struct semid_ds *buf; // параметры массива семафоров в целом
ushort *array;        // массив значений семафоров
}

*/

int main(int argc, char **argv)
{		
	int key, i, pid, p,  status, semid, shmid;
	int *shmaddr; // указатель на буфер с shared memory
	struct sembuf sops;
	p=atoi(argv[1]);             // кол-во сын. процессов
	key=ftok("FILE",'S');        //получаем уникальный ключ	
// создаем один семафор с определенными правами доступа
	semid=semget(key,1,0666|IPC_CREAT);
// создаем разделяемую память размером 1 элемент типа int
	shmid=shmget(key,1,0666|IPC_CREAT); 
// подключаемся к памяти
	shmaddr=shmat(shmid,NULL,0);
// инициализируем семафор значением 0
	semctl(semid,0, SETVAL, (int ) 0);
	sops.sem_num=0;  // чтобы использовать атомарную операцию
	sops.sem_flg=0;  // semop(semid,&sops,1)
//	cnt=*shmaddr;
	*shmaddr=0;      // обнуляем
	for (i=0; i<p; i++) {
		if((pid=fork())==-1) {
			perror("fork");
			exit(1);
		}
		if(pid==0) {
			// подключаемся к памяти
        		shmaddr=shmat(shmid,NULL,0);
//			key=ftok("FILE",'S');    // 
//        		semid=semget(key,1,0666|IPC_CREAT);
			sops.sem_num=0;
		        sops.sem_flg=0; 
				
			sops.sem_op=-i; 	// ожидание на семафоре
			semop(semid,&sops,1);   // i-го сына
			*shmaddr=*shmaddr+i*10;
			printf("%u процесс значение %d\n", i,*shmaddr);
			shmdt(shmaddr) ; // отключаемся от разделяемой памяти
//			sleep(i); //тест, что отец ждет обнуления семафора
			exit (0);
		}
	}
	i=0;
	while((wait(&status))!=-1) {
		sops.sem_op = 0;
                semop(semid, &sops, 1); //ожидание обнуления семафора 
		i++;
		sops.sem_op=i;      
                semop(semid,&sops,1);  // 
//		sops.sem_op = 0; 	
//		semop(semid, &sops, 1); //ожидание обнуления семафора 
	}
	shmdt(shmaddr) ; // отключаемся от разделяемой памяти
	shmctl(shmid, IPC_RMID, NULL); // уничтожаем разделяемую память
	semctl(semid,0,IPC_RMID, (int) 0);  //уничтожение семафора
	printf("Отец %d\n",0);
	return 0;
}
