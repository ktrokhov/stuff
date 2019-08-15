/* Программе в аргументах командной строки задаются: FILE - имя входного
бинарного файла и Р -кол-во процессов для обработки данных. Входной 
бинарный файл содержит 32-битные целые беззнаковые числа в представлении
Little-Endian (LE). Главный процесс должен создать Р процессов, которые
должны параллельно обработать файл FILE. Процесс с номером i (считая, что
процессы нумеруются от 0) вычисляет сумму чисел по модулю 2^32 на позициях
i, P+i, 2*P+i... в файле (позиции нумеруются от 0). После окончания 
обработки файла все процессы выводят на stdout по очереди в порядке 
возрастания номеров процессов накопившуюся у каждого процесса сумму.
Родитель дожидается завершения всех созданных им процессов, выводит на 
sddout число 0 и сам завершает свою работу с кодом завершения 0. Для 
синхронизации процессов использовать семафоры. */

/*LE - Порядок от младшего к старшему (т.е. в 0 бите младший разряд числа, в 
процессорах Intel X86 используется LE; 
  Запуск проги: ./com4 FILE 3 , например. FILE записан с помощью ./bineasy */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/sem.h>


/* ORIGIN       ОПРЕДЕЛЕНИЕ     
  SEEK_SET     начало файла.
  SEEK_CUR   текущая позиция указателя на файл.
  SEEK_END     конец файла */



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

union semun {
    int val;
    struct semid_ds *buf;
    ushort *array;
};
*/

int main(int argc, char **argv)  // Для >= 4 процессов
{		
	int key, i, pid, p, fd, ln, status, semid;
	struct sembuf sops;
	fd=open(argv[1], O_RDONLY ); // открыли бинарный файл 
	ln=lseek(fd, 0L, SEEK_END ); //длина файла в байтах
	p=atoi(argv[2]);             // кол-во сын. процессов
	key=ftok("FILE",'S');        //получаем уникальный ключ	
// создаем один семафор с определенными правами доступа
	semid=semget(key,1,0666|IPC_CREAT);
// инициализируем семафор значением 0
	semctl(semid,0, SETVAL, (int ) 0);
	sops.sem_num=0;  // чтобы использовать атомарную операцию
	sops.sem_flg=0;  // semop(semid,&sops,1)
	for (i=0; i<p; i++) {
		if((pid=fork())==-1) {
			perror("fork");
			exit(1);
		}
		if(pid==0) {
			int j=0 ,cur=0, tmp=0;
			for (j=0; ((i+j*p)<(ln/4)); j++) { // - ln/4 кол-во слов	
				lseek(fd,(i+j*p)*sizeof(int),SEEK_SET);
				read(fd,&cur,sizeof(int));
				tmp=tmp+cur;
				tmp=tmp&4294967295; // сумма по модулю 2^32
			}
			close(fd);
			cur=tmp;
//			key=ftok("FILE",'S');    // 
//        		semid=semget(key,1,0666|IPC_CREAT);
			sops.sem_num=0;
		        sops.sem_flg=0; 
				
			sops.sem_op=-i; 	// ожидание на семафоре
			semop(semid,&sops,1);   // i-го сына
			printf("%u процесс сумма %u\n", i,cur);
			sleep(i); //тест, что отец ждет обнуления семафора
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
	semctl(semid,0,IPC_RMID, (int) 0);  //уничтожение семафора
	close(fd);
	printf("Отец %d\n",0);
	return 0;
}
