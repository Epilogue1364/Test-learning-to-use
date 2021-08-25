/*如图所示，有一横截面为L1*L2的无限长度钢材，其四边边界温度分别为T1，T2，T3，T4
钢材密度ρ=7860kg/m3，导热系数k=40W/(m·K)，恒压热容Cp=460J/(kg·K)，无源项
求解该钢材温度分布可将其视为一个二维平面进行分析
假设L1=1mm，L2=1mm；T1=50℃；T2=50+100sin(pi*x/L2)；T3=50℃；T4=20℃。*/
#include <stdio.h>
#include<math.h>
#define pi 3.1415926535
#define ROW 100 //最大行
#define COL 100 //最大列
void create_Ab(double[][COL], double[], int, int);
void show_matrix(double[][COL], int, int);
void guass_elimination(double *[ROW], int, int);
void exchange_row(double *[ROW], int, int, int);
void show_solution(double *[ROW], int, int);
void back_substitution(double *[ROW], int, int);
void back_substitution(double *[ROW], int, int, int);

double x[COL];  // 存储唯一解x

int main()
{
	//矩阵计算变量
	double Receptacle[ROW][COL];    // 储存矩阵A
	double Vector[ROW];             // 储存向量b
	double *Ab_pointer[ROW];        // 储存增广矩阵(A,b)
	//double x[COL];
	int row, col;                   // 储存矩阵的行和列;
	int i;

	//数值分析变量
	double L1, L2, T1, T2, T3, T4, t, p, k, cp, dx, dy, tmax, dt;
	double ap0, aw, ae, as, an, ap, at, At, Ap;
	int Nx, Ny, Nt, NT;
	p = 7860;//密度
	k = 40;//导热系数
	cp = 460;//恒压热容
	T1 = 50;//左边界温度
	T3 = 50;//右边界温度
	T4 = 20;//下边界温度
	L1 = 0.001;//x方向 宽度
	L2 = 0.001;//y方向 高度
	dx = 0.0002;//delta x   为更快速地验证计算结果，减小计算量
	dy = 0.0002;//delta y
	Nx = L1 / dx + 1;//横向节点数 用于定义数组大小
	Ny = L2 / dy + 1;//纵向节点数
	printf("横向节点数Nx=%d\n", Nx);
	printf("纵向节点数Ny=%d\n", Ny);
	printf("二维截面宽L1=%0.3fm\n", L1);
	printf("二维截面高L2=%0.3fm\n", L2);
	printf("dx=%0.4f\n", dx);
	printf("dy=%0.4f\n", dy);

	double t0[6][6], t1[6][6];//本想将数组定义为变量数组，但是c++不能定义变量数组，所以只能提前计算好节点数，再定义数组大小
	double zx[6], zy[6];//点对应xy坐标
	double b0[16], b1[16];//全隐格式温度求解-----Vector[ROW]
	double A[16][16];//系数矩阵-----改用Receptacle[ROW][COL]xxxx———还是要用A矩阵，再在循环中，没经历一次循环用A给Receptacle做个代换
	NT = (Nx - 2)*(Ny - 2);
	row = col = NT;
	zx[0] = 0;//原点为0
	zy[0] = 0;
	for (int m = 1; m<Nx; m++)//点对应x坐标，用于上边界温度定义
	{
		zx[m] = zx[m - 1] + dx;
	}
	for (int n = 1; n<Ny; n++)//点对应y坐标，本次仿真中没有用到该值
	{
		zy[n] = zy[n - 1] + dy;
	}

	//定义全部节点初始值
	//从左至右，从上到下初始化内部   m=0表示上边界 n=0表示左边界

	for (int m = 0; m<Ny; m++)
	{
		for (int n = 0; n<Nx; n++)
		{
			t0[m][n] = T4;//钢材温度初始值同下边界温度，为20℃
		}
	}
	//初始上下边界条件
	for (int n = 0; n<Nx; n++)//初始化上下边
	{
		t0[0][n] = 50 + 100 * sin(pi*zx[n] / L1);//上边界
		t0[Ny - 1][n] = T4;//下边界
		t1[0][n] = 50 + 100 * sin(pi*zx[n] / L1);//未来值的边界不变
		t1[Ny - 1][n] = T4;//未来值的边界不变
	}
	//初始左右边界条件
	for (int m = 1; m<Ny; m++)//初始化左右边
	{
		t0[m][0] = T1;//左边
		t0[m][Nx - 1] = T3;//右边
		t1[m][0] = T1;//左边
		t1[m][Nx - 1] = T3;//右边
	}
	printf("输出初始值\n");
	for (int m = 0; m<Ny; m++)//输出初始值，检查有无错误
	{
		for (int n = 0; n<Nx; n++)
			printf("%10f ", t0[m][n]);
		printf("\n");
	}

	t = 0.02;//提前给定t，为加快测试
	printf("计算时间t=%0.4fs\n", t);
	at = k;//本网格划分中dx，dy大小相等，且每步积分四周的四个点到p点的距离与dx，dy大小相等，所以相互约去，只剩下导热系数k
	aw = ae = as = an = at;
	dt = 0.0005;//为加快测试进度，提前给定合适的dt
	ap0 = p*cp*dx*dy / dt;
	ap = ap0 + aw + ae + as + aw;//隐式格式中的ap表达式
	At = -at / ap0;//系数矩阵中的四周点系数
	Ap = ap / ap0;//系数矩阵中的求解点系数
	//ap = ap0;//显式格式中的ap与ap0大小相等
	Nt = t / dt;//迭代次数
	printf("步长dt=%0.4f\n", dt);
	printf("迭代次数Nt=%d\n", Nt);
	printf("ap0=%0.8f\n", ap0);
	printf("ap=%0.8f\n", ap);
	printf("Ap=%0.8f\n", Ap);
	printf("At=%0.8f\n", At);
	printf("隐式格式计算二维非稳态传热\n");
	//将求解域的初始值 也就是 (2,2)-(5,5)， 对应mn(1,1)-(4,4)
	printf("输出求解区域初始列矩阵\n");
	int r = 0;
	for (int m = 1; m<Ny - 1; m++)//输出初始值，检查有无错误
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
	///////////////////输出初始值，检查有无错误
	//////Ax=b中的b值获得
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
			printf("%10f ", Vector[r]);//Ax=b中的b值
			r++;
		}
	}
	r = 0;
	printf("\n");

	for (int m = 0; m<NT; m++)//输出初始值，检查有无错误
	{
		printf("%8.3f", x[m]);
	}
	printf("\n");
	
	for (int m = 0; m<NT; m++)//系数矩阵初始置零
	{
		for (int n = 0; n < NT; n++)
			A[m][n] = 0;
	}

	for (int m = 0; m<NT; m++)//系数矩阵赋值
	{
		for (int n = 0; n < NT; n++)
		{
			if (m == n)
				A[m][n] = Ap;//需替换成Ap，先用1检测
			if ((n == (m + 4)) || (n == (m - 4))|| (n == (m - 1))|| (n == (m + 1)))
				A[m][n] = At;
			if (((m == 3) && (n == 4)) || ((m == 4) && (n == 3)) || ((m == 7) && (n == 8)) || ((m == 8) && (n == 7)) || ((m == 11) && (n == 12)) || ((m == 12) && (n == 11)))
				A[m][n] = 0;
			printf("%0.8f ", A[m][n]);
			if (n == 3 || n == 7 || n == 11)
				printf("|");//Ax=b中的b值
		}
		printf("\n");
		if (m == 3 || m == 7 || m == 11 || m == 15)
			printf("-------------------------------------------------\n");//Ax=b中的b值
	}

	/////////////////////////迭代计算//////////////////////////////

	for (int u = 0; u < Nt; u++)
	{
		printf("输出b值\n");
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
				printf("%10f ", Vector[r]);//Ax=b中的b值
				r++;
			}
		}
		r = 0;
		for (int m = 0; m<NT; m++)//系数矩阵赋值
		{
			for (int n = 0; n < NT; n++)
			{
				Receptacle[m][n] = A[m][n];
				/*
				printf("%0.8f ", Receptacle[m][n]);
				if (n == 3 || n == 7 || n == 11)
					printf("|");//Ax=b中的b值
				*/
			}
			/*
			printf("\n");
			if (m == 3 || m == 7 || m == 11 || m == 15)
				printf("-------------------------------------------------\n");//Ax=b中的b值
			*/
		}
		create_Ab(Receptacle, Vector, row, col);
		printf("\n线性方程组的增广矩阵形式为:\n");
		//show_matrix(Receptacle, row, col + 1);  //输出增广矩阵（A，b）
		for (i = 0; i < ROW; i++)               // Ab_pointer储存增广矩阵(A,b)的每一行  make every Ab_pointer points every row of the augment matrix (A,b)
			Ab_pointer[i] = Receptacle[i];
		guass_elimination(Ab_pointer, row, col + 1);    // 用高斯消去法得到结果
		printf("检测返回x的值\n");
		for (int m = 0; m<NT; m++)//输出初始值，检查有无错误
		{
			printf("%8.3f", x[m]);
			if (m == 3 || m == 7 || m == 11 || m == 15)
				printf("\n-------------------------------------------------\n");//Ax=b中的b值
		}
		printf("\n");
		//printf("\n输入系数矩阵的大小(小于%d * %d, q退出): ", ROW, COL - 1);
	}
	for (int m = 1; m<Ny - 1; m++)
	{
		for (int n = 1; n < Nx - 1; n++)
		{
			t0[m][n] = x[r];
			
			printf("%10f ", t0[m][n]);//Ax=b中的b值
			r++;
		}
	}
	r = 0;
	printf("\n计算结果如下，其中第一行和第一列分别表示横坐标和纵坐标\n");
	printf("%10f ", zx[0]);

	for (int n = 0; n < Nx; n++)//输出x坐标
	{
		printf("%10f ", zx[n]);
	}
	printf("\n");
	for (int m = 0; m<Ny; m++)
	{
		printf("%10f ", zy[m]);//第一个列为列坐标
		for (int n = 0; n < Nx; n++)
		{
			printf("%10f ", t0[m][n]);//输出t时间的温度分布
		}
		printf("\n");
	}
	scanf("%lf", &L1);//防止控制台关闭
}


//子函数

/*************************************************
Function:       create_Ab
Description:    创造增广矩阵
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
Description:    输出矩阵
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
Description:    将矩阵变换为上三角矩阵，并通过高斯消去得到解!
*************************************************/
void guass_elimination(double *matrix[ROW], int row, int col)
{
	int result;         // 解决方案编号的标志
	int rankA, rankAb;  // 存储矩阵A或增广矩阵(A,b)的秩
	int i, j, k, l;
	double coe;         // 临时存储两行的倍数

						// 将矩阵变换为上三角矩阵
	for (i = 0; i < row; i++)
	{
		// 求交换矩阵行后的第一个非零值
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

		// 进行消除
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

	// 获得rank(A)
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

	// 获得rank(A,b)
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

	// rank(A)!=rank(A,b) => 没有解
	if (rankA != rankAb)
	{
		result = -1;
		printf("\n消除后:\n");
		show_solution(matrix, row, col);
		printf("\n线性方程是没有解的!\n");
	}

	// rank(A)=rank(A,b)=col => 唯一解
	else if (rankA == col - 1)
	{
		result = 0;
		printf("\n消除后:\n");
		//show_solution(matrix, row, col);
		printf("\n线性方程组只有一个解!\n");
		back_substitution(matrix, row, col - 1);
	}

	// rank(A)=rank(A,b)<col => 无穷解
	else
	{
		result = 1;
		printf("\n消除后:\n");
		show_solution(matrix, row, col);
		printf("\n线性方程组有无穷个解!\n");
		back_substitution(matrix, row, col - 1, rankA);
	}


	return;
}

/*************************************************
Function:       exchange_row
Description:    交换矩阵的行确保矩阵x[i][j]!=0
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
Description:    消除后输出上三角矩阵
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
Description:    通过反向代换找到唯一的解
*************************************************/
void back_substitution(double *matrix[ROW], int row, int col)
{
	int i, j;
	double temp;
	//double x[COL];  // 存储唯一解x

	for (i = 0; i < col; i++)
	{
		temp = matrix[col - 1 - i][col];
		for (j = 0; j < i; j++)
		{
			temp -= x[col - 1 - j] * matrix[col - 1 - i][col - 1 - j];
		}
		x[col - 1 - i] = temp / matrix[col - 1 - i][col - 1 - i];
	}

	// 输出结果
	printf("结果是:[", i, x[i]);
	for (i = 0; i < col; i++)
	{
		printf("%8.3f", x[i]);
	}
	printf("]\n");

	return;
}

/*************************************************
Function:       back_substitution
Description:    求无穷解
				第一步:求Ax=b的特殊解Xp
				第二步:求Ax=0的一般解Xn
*************************************************/
void back_substitution(double *matrix[ROW], int row, int col, int rankA)
{
	int i, j, k, n;
	int pivot[COL], free[COL];   //存储枢轴/自由位置
	double temp;
	double x_p[COL];            // 存储Ax=b的特殊解决方案Xp
	double x_n[COL][COL] = { 0 };  // 存储Ax=0的一般解

								   // 找出枢轴的位置
	for (i = 0; i < rankA; i++)
		for (j = i; j<col; j++)
			if (fabs(matrix[i][j])>0.00001)
			{
				pivot[i] = j;
				break;
			}

	// 找到自由的位置
	j = n = 0;
	for (i = 0; i < col; i++)
		if (i == pivot[j])
		{
			j++;
		}
		else
		{
			free[n] = i;
			x_p[i] = 1;   // 设置Xp的自由值
			n++;
		}

	// 求Ax=b的一个特解
	for (i = 0; i < rankA; i++)
	{
		n = rankA - 1 - i;    // 获取当前行号
		temp = matrix[n][col];
		for (j = pivot[n] + 1; j < col; j++)
		{
			temp -= x_p[j] * matrix[n][j];
		}
		x_p[pivot[n]] = temp / matrix[n][pivot[n]];
	}

	// 设置Xn的自由值
	for (i = 0; i < col - rankA; i++)
		x_n[i][free[i]] = 1;

	// 求出Ax=0的通解
	for (k = 0; k < col - rankA; k++)
	{
		for (i = 0; i < rankA; i++)
		{
			n = rankA - 1 - i;    // 获取当前行号
			temp = 0;
			for (j = pivot[n] + 1; j < col; j++)
			{
				temp -= x_n[k][j] * matrix[n][j];
			}
			x_n[k][pivot[n]] = temp / matrix[n][pivot[n]];
		}
	}

	// 输出结果
	printf("结果是 X=Xp+Xn.\n");
	printf("向量是 Xp is:[");
	for (i = 0; i < col; i++)
	{
		printf("%8.3f", x_p[i]);
	}
	printf("]\n");

	for (k = 0; k < col - rankA; k++)
	{
		printf("向量 Xn%d is:[", k + 1);
		for (i = 0; i < col; i++)
		{
			printf("%8.3f", x_n[k][i]);
		}

		printf("]\n");
	}

	return;
}