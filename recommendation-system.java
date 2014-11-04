package assignment;

import java.util.Random;

public class Recommendation {
		
	static Random random = new Random();
	
	static double[][] M = {	{0.0,0.0},
							{0.0,0.0}};
			
	public static void main(String[] args) {		
		System.out.println(recommend(M, 4)[0]);
	}
		
	public static String matrixToString(double[][] A) {
		String newLine = System.getProperty("line.separator");
		String s = "{" + newLine;
		int t1 = A.length;
		
		for(int i = 0; i<t1; ++i) {
			s = s + "{";
			int t2 = A[i].length;
			
			for(int j = 0; j < t2; ++j) {
				if(j<t2-1) {
					s = s + (double) A[i][j] + ",";
				} else {
					s = s + (double) A[i][j];
				}
			}
			
			if(i < t1-1) {
				s = s + "}," + newLine;
			} else {
				s = s + "}," + newLine;
			}
		}
		s = s + "};";	
		return s;
	}
	
	public static boolean isMatrix( double[][] A ) {
		if (A == null) {
			return false;
		} else {
			boolean notEmpty = isMatrixNotEmpty(A); 
			if(!notEmpty) {
				return false;
			} else {
				boolean rowsId = isMatrixRowsId(A);
				if(rowsId) {
					return true;
				}
				return false;
			}
						
		}	
	}
	
	public static boolean isMatrixNotEmpty(double[][] A) {
		int t = A.length;
		if (t == 0) {
			return false;
		}
		return true;
	}
	
	public static boolean isMatrixRowsId(double[][] A) {
		boolean e2 = true;
		
		int t = A.length;
		int t2 = A[0].length;
		for(int i = 1; i < t; ++i) {
			int a = A[i].length;
			if(a != t2) {
				e2 = false;
				break;
			} else {
				e2 = true;
			}
		}
		return e2;
	}
	
	public static double[][] multiplyMatrix(double[][] A, double[][] B) 
	{
		if(A[0].length!=B.length) {
			return null;
		} else {
			double[][] C = new double [A.length][B[0].length];
			
			for(int i = 0; i < A.length; ++i) {
				for(int j=0; j < B[0].length; ++j) {
					double pij = 0.0;
					for(int x=0; x<A[0].length; ++x) {
						pij = pij + (A[i][x])*(B[x][j]);
					}
					C[i][j] = pij;
				}
			}
			return C;
		}	
	}
	
	public static double[][] createMatrix( int n, int m, int k, int l) {
		double[][] R = new double [n][m];

		if(k > l || m == 0 || n == 0) {
			return null;
		}
		Random random = new Random();
			
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < m; ++j) {
				R[i][j] = (l-k) * random.nextDouble() + k;
			}
		}
		return R;
	}
	
	public static double rmse(double[][] M, double[][] P) {
		int counter = 0;
		double sum = 0;
		int mt = M.length;
		boolean matrix = checkMatrixDimensions(M,P);

		if(!matrix) {
			return -1;
		}

		for(int i = 0; i<mt; ++i) {
			int mtt = M[i].length;
			for(int j = 0; j<mtt; ++j) {
				if(M[i][j] == 0) {
					sum = sum +0;
				} else {
					double sqrd = (M[i][j] - P[i][j]);
					double sqrd2 = sqrd * sqrd;
					sum = sum + sqrd2;
					++counter;
				}		
			}
		}
		double average = sum/counter;
		double rmse = Math.sqrt(average);
		return rmse;
	}
	
	public static boolean checkMatrixDimensions(double M[][], double P[][]) {
		if(M.length!=P.length) {
			return false;
		}
		for(int i = 0; i < M.length; ++i) {
			if(M[i].length != P[i].length) {
				return false;
			}
		}
		return true;
	}
	
	public static double updateUElem( double[][] M, double[][] U, double[][] V, int r, int s ) {
		double Urs = 0;
		double sumV = 0, sumV2 = 0, sumUV = 0;
		int m = M[0].length; 
		int d = U[0].length;
		for(int j = 0; j < m; ++j) {
			if(M[r][j] != 0) {
				sumV2 = sumV2 + V[s][j]*V[s][j]; 
				sumUV = 0;
				for(int k = 0; k < d; ++k) {
					if(k != s) {
						sumUV = sumUV + U[r][k] * V[k][j];
					}
				}
				sumV = sumV + V[s][j] * (M[r][j]-sumUV);
			}
		}
		if(sumV2 != 0) {
			Urs = sumV/sumV2;
		} else {
			Urs = sumV;
		}
		return Urs;
	}
	
	public static double updateVElem( double[][] M, double[][] U, double[][] V, int r, int s ) {
		double Vrs = 0;
		double sumU = 0, sumU2 = 0, sumUV = 0;
		int n = M.length;
		int d = U[0].length;
		for(int i = 0; i < n; ++i) {
			if(M[i][s] != 0) {
				sumU2 = sumU2 + U[i][r]*U[i][r];
				sumUV = 0;
				for(int k=0; k<d; ++k) {
					if(k!=r) {
						sumUV = sumUV + U[i][k] * V[k][s];
					}
				}
				sumU = sumU + U[i][r] * (M[i][s] - sumUV);
			}
		}
		
		if(sumU2!=0) {
			Vrs = sumU/sumU2;		
		} else {
			Vrs = sumU;
		}
		return Vrs;		
	}
	
	public static double[][] initializeU (double[][] N, int d) {
		double sum = 0;
		double[][] U2 = new double [N.length][d];
		int t1 = N.length;
		int t2 = N[0].length;
		int n = 0;
		for(int i = 0; i < t1; ++i) {
			for(int j = 0; j<t2; ++j) {
				sum = sum + N[i][j];
				if(N[i][j]!=0) {
					++n;
				}
			}
		}
		int a = N[0].length;
		double average = sum/n;
		double v = Math.sqrt(average/a);
		for(int i = 0; i < t1; ++i) {
			for(int j = 0; j<d; ++j) {
				U2[i][j]=v+random.nextDouble();
			}
		}
		return U2;
	}
	
	public static double[][] initializeV (double[][] N, int d) {
		double sum = 0;
		double[][] V2 = new double [d][N[0].length];
		int t1 = N.length;
		int t2 = N[0].length;
		int n=0;
		for(int i = 0; i<t1; ++i) {
			for(int j = 0; j<t2; ++j) {
				sum = sum + N[i][j];
				if(N[i][j]!=0) {
					++n;
				}
			}
		}
		int a = N[0].length;
		double average = sum/n;
		double v = Math.sqrt(average/a);
		for(int i = 0; i<d; ++i) {
			for(int j = 0; j<t2; ++j) {
				V2[i][j] = v+random.nextDouble();
			}
		}
		return V2;
	}
	
	public static double[][] optimizeU( double[][] M, double[][] U, double[][] V) {
		double[][] newU = U;
		for(int i = 0; i < U.length; ++i)
		{
			for(int j = 0; j < U[0].length; ++j)
			{
				newU[i][j] = updateUElem(M, newU, V, i, j);
			}
		}
		return newU;	
	}

	public static double[][] optimizeV( double[][] M, double[][] U, double[][] V) {
		double[][] newV = V;
		for(int i = 0; i < V.length; ++i) {
			for(int j = 0; j < V[0].length; ++j) {
				newV[i][j] = updateVElem(M, U, newV, i, j);
			}
		}
		return newV;		
	}
	
	public static int[] recommend( double[][] M, int d) {
		int[] I = new int [M.length];
		double NewRMSE = 0;
		double[][] U = initializeU(M, d);
		double[][] V = initializeV(M, d);
		double[][] P = new double[U.length][V[0].length];
		double[][] Uintermediaire = new double[M.length][d];
		double[][] Vintermediaire = new double[d][M[0].length];
		P = multiplyMatrix(U, V);
		double rmse = rmse(M, P);
		int repetition = nbrRepetition(M);

		for(int i = 0; i < repetition; ++i) {
			Uintermediaire = initializeU(M, d);
			Vintermediaire = initializeV(M, d);
			P = multiplyMatrix(Uintermediaire, Vintermediaire);
			rmse = rmse(M, P);
			if(rmse<NewRMSE) {
				NewRMSE = rmse;
				U = Uintermediaire;
				V = Vintermediaire;
			}
		}
		
		double rmse1 = 0, rmse2 = 0;

		if(!isMatrix(M)) {
			for(int i=0; i<M.length; ++i) {
				I[i]=-1;
			}
		} else {
			do {
				rmse1=rmse2;
				U = optimizeU(M, U, V);
				V = optimizeV(M, U, V);
				P = multiplyMatrix(U, V);
				rmse2 = rmse(M, P);
			} while(Math.abs(rmse1-rmse2)>0.000001);
		}
		
		int R[][] = new int[M.length][M[0].length];

		for(int i = 0; i < M.length; ++i) {
			for(int j = 0; j < M[0].length; ++j) {
				if(M[i][j] == 0) {
					R[i][j] = 1;
				} else if(M[i][j]!=0) {
					R[i][j] = -1;
				}
			}
		}
		
		for(int i = 0; i < M.length; ++i) {
			boolean full = false;
			int firstVal = R[i][0];
			for(int j = 0; j < M[0].length; ++j) {
				if(R[i][j]! = firstVal) {
					full = true;
				}
			}
			if(full) {
				I[i] = 0;
			} else {
				I[i] = -1;
			}
		}
		
		for(int i = 0; i < P.length; ++i) {
			if(I[i] != -1) {
				double max = 0;
				for(int j = 0; j < P[0].length; ++j) {
					if(R[i][j] == 1) {
						if(P[i][j] > max) {
							max = P[i][j];
							I[i] = j;
						}
					}
				}
			}
		}
		return I;
	}
	
	public static int nbrRepetition(double M[][]) {
		int tLine = M.length;
		int tColumn = M[0].length;
		double matrix = (tLine + tColumn)/2;
		double log = (Math.log(matrix)/Math.log(Math.E));
		int factor = (int) (1000000/Math.pow(Math.E, log));
		return factor;
	}
}