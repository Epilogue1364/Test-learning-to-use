/*��ͼ��ʾ����һ�����ΪL1*L2�����޳��ȸֲģ����ı߽߱��¶ȷֱ�ΪT1��T2��T3��T4
�ֲ��ܶȦ�=7860kg/m3������ϵ��k=40W/(m��K)����ѹ����Cp=460J/(kg��K)����Դ��
���øֲ��¶ȷֲ��ɽ�����Ϊһ����άƽ����з���
����L1=1mm��L2=1mm��T1=50�棻T2=50+100sin(pi*x/L2)��T3=50�棻T4=20�档*/
#include <stdio.h>
#include<math.h>
#define pi 3.1415926535
#define ROW 100 //�����
#define COL 100 //�����
void create_Ab(double[][COL], double[], int, int);
void show_matrix(double[][COL], int, int);
void guass_elimination(double *[ROW], int, int);
void exchange_row(double *[ROW], int, int, int);
void show_solution(double *[ROW], int, int);
void back_substitution(double *[ROW], int, int);
void back_substitution(double *[ROW], int, int, int);

double x[COL];  // �洢Ψһ��x

int main()
{
	//����������
	double Receptacle[ROW][COL];    // �������A
	double Vector[ROW];             // ��������b
	double *Ab_pointer[ROW];        // �����������(A,b)
	//double x[COL];
	int row, col;                   // ���������к���;
	int i;

	//��ֵ��������
	double L1, L2, T1, T2, T3, T4, t, p, k, cp, dx, dy, tmax, dt;
	double ap0, aw, ae, as, an, ap, at, At, Ap;
	int Nx, Ny, Nt, NT;
	p = 7860;//�ܶ�
	k = 40;//����ϵ��
	cp = 460;//��ѹ����
	T1 = 50;//��߽��¶�
	T3 = 50;//�ұ߽��¶�
	T4 = 20;//�±߽��¶�
	L1 = 0.001;//x���� ���
	L2 = 0.001;//y���� �߶�
	dx = 0.0002;//delta x   Ϊ�����ٵ���֤����������С������
	dy = 0.0002;//delta y
	Nx = L1 / dx + 1;//����ڵ��� ���ڶ��������С
	Ny = L2 / dy + 1;//����ڵ���
	printf("����ڵ���Nx=%d\n", Nx);
	printf("����ڵ���Ny=%d\n", Ny);
	printf("��ά�����L1=%0.3fm\n", L1);
	printf("��ά�����L2=%0.3fm\n", L2);
	printf("dx=%0.4f\n", dx);
	printf("dy=%0.4f\n", dy);

	double t0[6][6], t1[6][6];//���뽫���鶨��Ϊ�������飬����c++���ܶ���������飬����ֻ����ǰ����ýڵ������ٶ��������С
	double zx[6], zy[6];//���Ӧxy����
	double b0[16], b1[16];//ȫ����ʽ�¶����-----Vector[ROW]
	double A[16][16];//ϵ������-----����Receptacle[ROW][COL]xxxx����������Ҫ��A��������ѭ���У�û����һ��ѭ����A��Receptacle��������
	NT = (Nx - 2)*(Ny - 2);
	row = col = NT;
	zx[0] = 0;//ԭ��Ϊ0
	zy[0] = 0;
	for (int m = 1; m<Nx; m++)//���Ӧx���꣬�����ϱ߽��¶ȶ���
	{
		zx[m] = zx[m - 1] + dx;
	}
	for (int n = 1; n<Ny; n++)//���Ӧy���꣬���η�����û���õ���ֵ
	{
		zy[n] = zy[n - 1] + dy;
	}

	//����ȫ���ڵ��ʼֵ
	//�������ң����ϵ��³�ʼ���ڲ�   m=0��ʾ�ϱ߽� n=0��ʾ��߽�

	for (int m = 0; m<Ny; m++)
	{
		for (int n = 0; n<Nx; n++)
		{
			t0[m][n] = T4;//�ֲ��¶ȳ�ʼֵͬ�±߽��¶ȣ�Ϊ20��
		}
	}
	//��ʼ���±߽�����
	for (int n = 0; n<Nx; n++)//��ʼ�����±�
	{
		t0[0][n] = 50 + 100 * sin(pi*zx[n] / L1);//�ϱ߽�
		t0[Ny - 1][n] = T4;//�±߽�
		t1[0][n] = 50 + 100 * sin(pi*zx[n] / L1);//δ��ֵ�ı߽粻��
		t1[Ny - 1][n] = T4;//δ��ֵ�ı߽粻��
	}
	//��ʼ���ұ߽�����
	for (int m = 1; m<Ny; m++)//��ʼ�����ұ�
	{
		t0[m][0] = T1;//���
		t0[m][Nx - 1] = T3;//�ұ�
		t1[m][0] = T1;//���
		t1[m][Nx - 1] = T3;//�ұ�
	}
	printf("�����ʼֵ\n");
	for (int m = 0; m<Ny; m++)//�����ʼֵ��������޴���
	{
		for (int n = 0; n<Nx; n++)
			printf("%10f ", t0[m][n]);
		printf("\n");
	}

	t = 0.02;//��ǰ����t��Ϊ�ӿ����
	printf("����ʱ��t=%0.4fs\n", t);
	at = k;//�����񻮷���dx��dy��С��ȣ���ÿ���������ܵ��ĸ��㵽p��ľ�����dx��dy��С��ȣ������໥Լȥ��ֻʣ�µ���ϵ��k
	aw = ae = as = an = at;
	dt = 0.0005;//Ϊ�ӿ���Խ��ȣ���ǰ�������ʵ�dt
	ap0 = p*cp*dx*dy / dt;
	ap = ap0 + aw + ae + as + aw;//��ʽ��ʽ�е�ap���ʽ
	At = -at / ap0;//ϵ�������е����ܵ�ϵ��
	Ap = ap / ap0;//ϵ�������е�����ϵ��
	//ap = ap0;//��ʽ��ʽ�е�ap��ap0��С���
	Nt = t / dt;//��������
	printf("����dt=%0.4f\n", dt);
	printf("��������Nt=%d\n", Nt);
	printf("ap0=%0.8f\n", ap0);
	printf("ap=%0.8f\n", ap);
	printf("Ap=%0.8f\n", Ap);
	printf("At=%0.8f\n", At);
	printf("��ʽ��ʽ�����ά����̬����\n");
	//�������ĳ�ʼֵ Ҳ���� (2,2)-(5,5)�� ��Ӧmn(1,1)-(4,4)
	printf("�����������ʼ�о���\n");
	int r = 0;
	for (int m = 1; m<Ny - 1; m++)//�����ʼֵ��������޴���
	{
		for (int n = 1; n < Nx - 1; n++)
		{
			x[r]= t0[m][n];
			printf("%8.3f", x[r]);
			r++;
		}
	}
	r = 0;
	printf("\n");
	///////////////////�����ʼֵ��������޴���
	//////Ax=b�е�bֵ���
	for (int m = 1; m<Ny-1; m++)
	{
		for (int n = 1; n < Nx - 1; n++) 
		{
			Vector[r] = x[r];
				if (m == 1)
					Vector[r] = Vector[r] + an*t0[m - 1][n] / ap0;
				if (m == 4)
					Vector[r] = Vector[r] + as*t0[m + 1][n] / ap0;
				if (n == 1)
					Vector[r] = Vector[r] + aw*t0[m][n - 1] / ap0;
				if (n == 4)
					Vector[r] = Vector[r] + ae*t0[m][n + 1] / ap0;
			printf("%10f ", Vector[r]);//Ax=b�е�bֵ
			r++;
		}
	}
	r = 0;
	printf("\n");

	for (int m = 0; m<NT; m++)//�����ʼֵ��������޴���
	{
		printf("%8.3f", x[m]);
	}
	printf("\n");
	
	for (int m = 0; m<NT; m++)//ϵ�������ʼ����
	{
		for (int n = 0; n < NT; n++)
			A[m][n] = 0;
	}

	for (int m = 0; m<NT; m++)//ϵ������ֵ
	{
		for (int n = 0; n < NT; n++)
		{
			if (m == n)
				A[m][n] = Ap;//���滻��Ap������1���
			if ((n == (m + 4)) || (n == (m - 4))|| (n == (m - 1))|| (n == (m + 1)))
				A[m][n] = At;
			if (((m == 3) && (n == 4)) || ((m == 4) && (n == 3)) || ((m == 7) && (n == 8)) || ((m == 8) && (n == 7)) || ((m == 11) && (n == 12)) || ((m == 12) && (n == 11)))
				A[m][n] = 0;
			printf("%0.8f ", A[m][n]);
			if (n == 3 || n == 7 || n == 11)
				printf("|");//Ax=b�е�bֵ
		}
		printf("\n");
		if (m == 3 || m == 7 || m == 11 || m == 15)
			printf("-------------------------------------------------\n");//Ax=b�е�bֵ
	}

	/////////////////////////��������//////////////////////////////

	for (int u = 0; u < Nt; u++)
	{
		printf("���bֵ\n");
		for (int m = 1; m<Ny - 1; m++)
		{
			for (int n = 1; n < Nx - 1; n++)
			{
				Vector[r] = x[r];
				if (m == 1)
					Vector[r] = Vector[r] + an*t0[m - 1][n] / ap0;
				if (m == 4)
					Vector[r] = Vector[r] + as*t0[m + 1][n] / ap0;
				if (n == 1)
					Vector[r] = Vector[r] + aw*t0[m][n - 1] / ap0;
				if (n == 4)
					Vector[r] = Vector[r] + ae*t0[m][n + 1] / ap0;
				printf("%10f ", Vector[r]);//Ax=b�е�bֵ
				r++;
			}
		}
		r = 0;
		for (int m = 0; m<NT; m++)//ϵ������ֵ
		{
			for (int n = 0; n < NT; n++)
			{
				Receptacle[m][n] = A[m][n];
				/*
				printf("%0.8f ", Receptacle[m][n]);
				if (n == 3 || n == 7 || n == 11)
					printf("|");//Ax=b�е�bֵ
				*/
			}
			/*
			printf("\n");
			if (m == 3 || m == 7 || m == 11 || m == 15)
				printf("-------------------------------------------------\n");//Ax=b�е�bֵ
			*/
		}
		create_Ab(Receptacle, Vector, row, col);
		printf("\n���Է���������������ʽΪ:\n");
		//show_matrix(Receptacle, row, col + 1);  //����������A��b��
		for (i = 0; i < ROW; i++)               // Ab_pointer�����������(A,b)��ÿһ��  make every Ab_pointer points every row of the augment matrix (A,b)
			Ab_pointer[i] = Receptacle[i];
		guass_elimination(Ab_pointer, row, col + 1);    // �ø�˹��ȥ���õ����
		printf("��ⷵ��x��ֵ\n");
		for (int m = 0; m<NT; m++)//�����ʼֵ��������޴���
		{
			printf("%8.3f", x[m]);
			if (m == 3 || m == 7 || m == 11 || m == 15)
				printf("\n-------------------------------------------------\n");//Ax=b�е�bֵ
		}
		printf("\n");
		//printf("\n����ϵ������Ĵ�С(С��%d * %d, q�˳�): ", ROW, COL - 1);
	}
	for (int m = 1; m<Ny - 1; m++)
	{
		for (int n = 1; n < Nx - 1; n++)
		{
			t0[m][n] = x[r];
			
			printf("%10f ", t0[m][n]);//Ax=b�е�bֵ
			r++;
		}
	}
	r = 0;
	printf("\n���������£����е�һ�к͵�һ�зֱ��ʾ�������������\n");
	printf("%10f ", zx[0]);

	for (int n = 0; n < Nx; n++)//���x����
	{
		printf("%10f ", zx[n]);
	}
	printf("\n");
	for (int m = 0; m<Ny; m++)
	{
		printf("%10f ", zy[m]);//��һ����Ϊ������
		for (int n = 0; n < Nx; n++)
		{
			printf("%10f ", t0[m][n]);//���tʱ����¶ȷֲ�
		}
		printf("\n");
	}
	scanf("%lf", &L1);//��ֹ����̨�ر�
}


//�Ӻ���

/*************************************************
Function:       create_Ab
Description:    �����������
*************************************************/
void create_Ab(double matrix[ROW][COL], double vector[ROW], int row, int col)
{
	int i;

	for (i = 0; i < row; i++)
		matrix[i][col] = vector[i];

	return;
}

/*************************************************
Function:       show_matrix
Description:    �������
*************************************************/
void show_matrix(double matrix[ROW][COL], int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
			printf("%-8.3f", matrix[i][j]);
		putchar('\n');
	}

	return;
}

/*************************************************
Function:       guass_elimination
Description:    ������任Ϊ�����Ǿ��󣬲�ͨ����˹��ȥ�õ���!
*************************************************/
void guass_elimination(double *matrix[ROW], int row, int col)
{
	int result;         // ���������ŵı�־
	int rankA, rankAb;  // �洢����A���������(A,b)����
	int i, j, k, l;
	double coe;         // ��ʱ�洢���еı���

						// ������任Ϊ�����Ǿ���
	for (i = 0; i < row; i++)
	{
		// �󽻻������к�ĵ�һ������ֵ
		for (j = i; j < col - 1; j++)
		{
			if (fabs(matrix[i][j])<0.00001)
			{
				exchange_row(matrix, i, j, row);
				if (fabs(matrix[i][j])>0.00001)
				{
					break;
				}
			}
			else
				break;
		}

		// ��������
		for (k = i + 1; k < row; k++)
		{
			if (matrix[i][j] == 0)
				break;
			coe = matrix[k][j] / matrix[i][j];
			for (l = j; l < col; l++)
			{
				matrix[k][l] -= coe*matrix[i][l];
			}
		}
	}

	rankA = rankAb = 0;

	// ���rank(A)
	for (i = 0; i < row; i++)
	{
		for (j = i; j<col - 1; j++)
		{
			if (fabs(matrix[i][j])>0.00001)
			{
				rankA++;
				break;
			}
		}
	}

	// ���rank(A,b)
	for (i = 0; i < row; i++)
	{
		for (j = i; j<col; j++)
		{
			if (fabs(matrix[i][j])>0.00001)
			{
				rankAb++;
				break;
			}
		}
	}

	// rank(A)!=rank(A,b) => û�н�
	if (rankA != rankAb)
	{
		result = -1;
		printf("\n������:\n");
		show_solution(matrix, row, col);
		printf("\n���Է�����û�н��!\n");
	}

	// rank(A)=rank(A,b)=col => Ψһ��
	else if (rankA == col - 1)
	{
		result = 0;
		printf("\n������:\n");
		//show_solution(matrix, row, col);
		printf("\n���Է�����ֻ��һ����!\n");
		back_substitution(matrix, row, col - 1);
	}

	// rank(A)=rank(A,b)<col => �����
	else
	{
		result = 1;
		printf("\n������:\n");
		show_solution(matrix, row, col);
		printf("\n���Է��������������!\n");
		back_substitution(matrix, row, col - 1, rankA);
	}


	return;
}

/*************************************************
Function:       exchange_row
Description:    �����������ȷ������x[i][j]!=0
(except matrix[i][j] ~ matrix[row-1][j] all 0)
*************************************************/
void exchange_row(double *matrix[ROW], int i, int j, int row)
{
	int k;
	double *temp;

	for (k = i + 1; k<row; k++)
	{
		if (fabs(matrix[k][j])>0.00001)
		{
			temp = matrix[i];
			matrix[i] = matrix[k];
			matrix[k] = temp;
			return;
		}
	}

	return;
}

/*************************************************
Function:       show_solution
Description:    ��������������Ǿ���
*************************************************/
void show_solution(double *matrix[ROW], int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
			printf("%-8.3f", matrix[i][j]);
		putchar('\n');
	}

	return;
}

/*************************************************
Function:       back_substitution
Description:    ͨ����������ҵ�Ψһ�Ľ�
*************************************************/
void back_substitution(double *matrix[ROW], int row, int col)
{
	int i, j;
	double temp;
	//double x[COL];  // �洢Ψһ��x

	for (i = 0; i < col; i++)
	{
		temp = matrix[col - 1 - i][col];
		for (j = 0; j < i; j++)
		{
			temp -= x[col - 1 - j] * matrix[col - 1 - i][col - 1 - j];
		}
		x[col - 1 - i] = temp / matrix[col - 1 - i][col - 1 - i];
	}

	// ������
	printf("�����:[", i, x[i]);
	for (i = 0; i < col; i++)
	{
		printf("%8.3f", x[i]);
	}
	printf("]\n");

	return;
}

/*************************************************
Function:       back_substitution
Description:    �������
				��һ��:��Ax=b�������Xp
				�ڶ���:��Ax=0��һ���Xn
*************************************************/
void back_substitution(double *matrix[ROW], int row, int col, int rankA)
{
	int i, j, k, n;
	int pivot[COL], free[COL];   //�洢����/����λ��
	double temp;
	double x_p[COL];            // �洢Ax=b������������Xp
	double x_n[COL][COL] = { 0 };  // �洢Ax=0��һ���

								   // �ҳ������λ��
	for (i = 0; i < rankA; i++)
		for (j = i; j<col; j++)
			if (fabs(matrix[i][j])>0.00001)
			{
				pivot[i] = j;
				break;
			}

	// �ҵ����ɵ�λ��
	j = n = 0;
	for (i = 0; i < col; i++)
		if (i == pivot[j])
		{
			j++;
		}
		else
		{
			free[n] = i;
			x_p[i] = 1;   // ����Xp������ֵ
			n++;
		}

	// ��Ax=b��һ���ؽ�
	for (i = 0; i < rankA; i++)
	{
		n = rankA - 1 - i;    // ��ȡ��ǰ�к�
		temp = matrix[n][col];
		for (j = pivot[n] + 1; j < col; j++)
		{
			temp -= x_p[j] * matrix[n][j];
		}
		x_p[pivot[n]] = temp / matrix[n][pivot[n]];
	}

	// ����Xn������ֵ
	for (i = 0; i < col - rankA; i++)
		x_n[i][free[i]] = 1;

	// ���Ax=0��ͨ��
	for (k = 0; k < col - rankA; k++)
	{
		for (i = 0; i < rankA; i++)
		{
			n = rankA - 1 - i;    // ��ȡ��ǰ�к�
			temp = 0;
			for (j = pivot[n] + 1; j < col; j++)
			{
				temp -= x_n[k][j] * matrix[n][j];
			}
			x_n[k][pivot[n]] = temp / matrix[n][pivot[n]];
		}
	}

	// ������
	printf("����� X=Xp+Xn.\n");
	printf("������ Xp is:[");
	for (i = 0; i < col; i++)
	{
		printf("%8.3f", x_p[i]);
	}
	printf("]\n");

	for (k = 0; k < col - rankA; k++)
	{
		printf("���� Xn%d is:[", k + 1);
		for (i = 0; i < col; i++)
		{
			printf("%8.3f", x_n[k][i]);
		}

		printf("]\n");
	}

	return;
}