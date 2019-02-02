#include<iostream>
#include<cstdlib>
#include<vector>
#include<cmath>
using namespace std;

template<class T>class MyMatrix;
template<class T> std::istream& operator>>(std::istream& in, MyMatrix<T>& dst);
template<class T> std::ostream& operator<<(std::ostream& out, const MyMatrix<T>& src);

template <class T>
class MyMatrix{
	const int row; //����
	const int column; //���������ھ������ɵ������Ĺ������ǲ���ģ������Ҫת�ã������½�һ������ 
	T* ptr; //��һ�����鱣����������Ԫ��
public:
	//���췽ʽ1������������������ָ�� 
	MyMatrix(int r, int c, T* p): row(r), column(c)
	{
		ptr = new T[row * column];
		int i;
		for(i = 0; i < row * column; i++)
		{
			ptr[i] = p[i];
		}
	}
	//���췽ʽ2�� ��������������������Ԫ����ʱ��Ϊ0,�����������������г�ʼ�� 
	MyMatrix(int r, int c): row(r), column(c)
	{
		ptr = new T[row * column];
		for(int i = 0; i < row * column; i++)
		{
			ptr[i] = 0;
		}
	}
	
	MyMatrix<T>(const MyMatrix<T>& A): row(A.row), column(A.column)
	{
		ptr = new T[row * column];
		int i;
		for(i = 0; i < row * column; i++)
		{
			ptr[i] = A.ptr[i];
		}
	}
	
	MyMatrix<T>& operator= (const MyMatrix<T>& rv)
	{
		if(this == &rv)
		{
			return *this;
		}
		else if((row != rv.row)||(column != rv.column))
		{
			std::cout<<"matrix set error: not same size!\n";
			return *this;
		}
		else
		{
			delete[] this->ptr;
			this->ptr = new T[row * column];
			int i;
			for(i = 0; i < row * column; i++)
			{
				this->ptr[i] = rv.ptr[i];
			}
		}
	}
	
	MyMatrix<T> operator+ (const MyMatrix<T>& rv) const
	{
		int i, j;
		if((row != rv.row)||(column != rv.column))
		{
			MyMatrix<T> A(1,1);
			return A;
		}
		else
		{
			MyMatrix<T> A(row, column);
			for(i = 0; i < row; i++)
			{
				 for(j = 0; j < column; j++)
				 {
				 	A.SetNumber(i, j, GetNumber(i, j) + rv.GetNumber(i, j));
				 }
			}
			return A;
		}
	}
	
	MyMatrix<T> operator- (const MyMatrix<T>& rv) const
	{
		int i, j;
		if((row != rv.row)||(column != rv.column))
		{
			MyMatrix<T> A(1,1);
			return A;
		}
		else
		{
			MyMatrix<T> A(row, column);
			for(i = 0; i < row; i++)
			{
				 for(j = 0; j < column; j++)
				 {
				 	A.SetNumber(i, j, GetNumber(i, j) - rv.GetNumber(i, j));
				 }
			}
			return A;
		}
	}
	
	//�ھ�������ɺ󽫵�i�е�j����Ϊnumber��������Ϊi,j����0��ʼ 
	void SetNumber(int i, int j, T number)
	{
		int position = i * column + j;
		if( i < 0 || j < 0 || i >= row || j >= column)
		{
			std::cout<<"SetNumber error: the element does not exist!\n";
			return;
		}
		else
		{
			ptr[position] = number;
			return;
		}
	}
	
	//��þ����i�С���j�е�Ԫ�أ�i,j��0��ʼ
	/*Ҫ���������������Ա��������GetNumber��SetNumber�򽻵�������һ������Ԫ�صĴ洢��ʽ
	 �����ı䣨����ϡ����󣩣�ֻ�ø��������ײ㺯������*/ 
	T GetNumber(int i, int j) const
	{
		int position = i * column + j;
		if(i < 0 || j < 0 || i >= row || j>= column )
		{
			std::cout<<"GetNumber error: element does not exist!\n";
			return 0;
		}
		else
		{
			return ptr[position];
		}
	}
	T infty_norm()
	{
		T max = (T)0;
		for(int i = 0; i < row; i++)
		{
			T temp = (T)0;
			for(int j = 0; j < column; j++)
			{
				T element = MyMatrix<T>::GetNumber(i, j);
				if(element >= (T)0)
				{
					temp += element;
				}
				else
				{
					temp -= element;
				}
			}
			if(temp > max)
			{
				max = temp;
			}
		}
		return max;
	}
	
	MyMatrix<T> operator* (const MyMatrix<T>& rv) const
	{
		if(column != rv.row)
		{
			MyMatrix<T> A(1,1);
			return A;
		}
		else
		{
			int i, j, k;
			T number;
			MyMatrix<T> A(row, rv.column);
			for(i = 0; i < row; i++)
			{
				for(j = 0; j < rv.column; j++)
				{
					number = 0;
					for(k = 0; k < column; k++)
					{
						number = number + GetNumber(i, k) * rv.GetNumber(k, j);
					}
					A.SetNumber(i, j, number);
				}
			}
			return A;
		}
	}
	friend std::istream& operator>> <>(std::istream& in, MyMatrix<T>& dst);
	friend std::ostream& operator<< <>(std::ostream& out, const MyMatrix<T>& src);
	//��һ�ֳ��ȱ任������i�г�k 
	void RowMutiply(int i, T k);
	//�ڶ��ֳ��ȱ任������j�г�k�ӵ���i��
	void RowAdd(int i, int j, T k); 
	//�����ֳ��ȱ任�� ����i�����j�л���
	void RowChange(int i, int j);
	//��˹��Ԫ������������ݻ��� ����rowchange�Ĵ��� 
	int GauessElimination(MyMatrix<T>& Mirror); 
	int GauessElimination();
	//�������� ���ڰ��е�˳������������Ԫ�� 
	class iterator;
	friend class iterator;
	class iterator{
		MyMatrix& s;
		int index;
	public:
		iterator(MyMatrix& st) : s(st), index(0) {}
		//to create the "end sentinel" iterator
		iterator(MyMatrix& st, bool) : s(st), index(s.row * s.column - 1){}
		T operator*() const
		{
			return s.ptr[index];
		}
		T operator++() //prefix form
		{
			if(index + 1 >= s.row * s.column)
			{
				std::cout<<"prefix ++ error: iterator move out of range!\n";
				return (T)0;
			}
			else
			{
				return s.ptr[++index];
			}
		}
		T operator++(int) //postfix form
		{
			if(index + 1 >= s.row * s.column)
			{
				std::cout<<"postfix ++ error: iterator move out of range!\n";
				return (T)0;
			}
			else
			{
				return s.ptr[index++];
			}
		}
		iterator& operator+=(int amount)
		{
			if(index + amount >= s.row * s.column)
			{
				std::cout<<"+= error: iterator trying to move out of the bound!\n";
				return *this;
			}
			else
			{
				index += amount;
				return *this;
			}
		}
		T operator--() //���ھ���������ԣ�������Ҫһ���������ƶ�һ��Ĳ�����
		{
			if(index + s.column >= s.row * s.column)
			{
				std::cout<<"prefix -- error: iterator cannot move down!\n";
				return (T)0;
			}
			else
			{
				index += s.column;
				return s.ptr[index];
			}
		}
		T operator--(int)
		{
			if(index + s.column >= s.row * s.column)
			{
				std::cout<<"prefix -- error: iterator cannot move down!\n";
				return (T)0;
			}
			else
			{
				index += s.column;
				return s.ptr[index - s.column];
			}
		}
		//to see if you are at the end
		bool operator==(const iterator& rv) const
		{
			return index == rv.index;
		}
		bool operator!=(const iterator& rv) const
		{
			return index != rv.index;
		}
		friend std::ostream& operator<<(std::ostream& os, const iterator& it)
		{
			return os << *it;
		}
	};
	iterator begin()
	{
		return iterator(*this);
	}
	//create the "end sential"
	iterator end()
	{
		return iterator(*this, true);
	}
	//�������� 
	virtual int GetRow() const
	{
		return row;
	}
	virtual ~MyMatrix()
	{
		delete[] ptr;
	}
};


template <class T> std::istream& operator>>(std::istream& in, MyMatrix<T>& dst)
{
	int i, j;
	T temp;
	for(i = 0; i < dst.row; i++)
	{
		for(j = 0; j < dst.column; j++)
		{
			in >> temp;
			dst.SetNumber(i, j, temp);
		}
	}
	return in;//Ҫ����in,��Ϊ����Ҫ����Ϊ�ö���ֵ 
}

template <class T> std::ostream& operator<<(std::ostream& out, const MyMatrix<T>& src)
{
	int i, j;
	for(i = 0; i < src.row; i++)
	{
		for(j = 0; j < src.column; j++)
		{
			out<<src.GetNumber(i, j)<<' ';
		}
		out<<'\n';
	}
	out<<'\n';
}

template<class T>
void MyMatrix<T>::RowMutiply(int i, T k)
{
	if(i < 0||i >= row)
	{
		std::cout<<"Row Multiply error: the row doesn't exist!\n";
		return;
	}
	else
	{
		int j;
		for(j = 0; j < column; j++)
		{
			SetNumber(i, j, k * GetNumber(i, j));
		}
		return;
	}
}

template<class T>
void MyMatrix<T>::RowAdd(int i, int j, T k)
{
	if(i < 0 || i >= row || j < 0 || j >= row)
	{
		std::cout<<"Row Add error: a row does not exist!\n";
		return;
	}
	else
	{
		int t;
		for(t = 0; t < column; t++)
		{
			SetNumber(i, t, GetNumber(i, t) + k * GetNumber(j, t));
		}
		return;
	}
}

template<class T>
void MyMatrix<T>::RowChange(int i, int j)
{
	if(i < 0 || i >= row || j < 0 || j >= row)
	{
		std::cout<<"Row Change error: a row does not exist!\n";
		return;
	}
	else
	{
		T temp;
		int t;
		for(t = 0; t < column; t++)
		{
			temp = GetNumber(i, t);
			SetNumber(i, t, GetNumber(j, t));
			SetNumber(j, t, temp);
		}
		return;
	}
}

template<class T>
T abs(T m)
{
	if(m >= (T)0)
	{
		return m;
	}
	else
	{
		return -m;
	}
}
template<class T>
int MyMatrix<T>::GauessElimination(MyMatrix<T>& Mirror)
{
	int result = 1;
	T temp;
	int MaxRow;
	T max;
	int i, j;
	int CurrentRow = 0;
	if(Mirror.GetRow() != row)
	{
		std::cout<<"Mirror matrix's size should be adjust!\n";
		return 0;
	}
	for( i = 0; i < column; i++)
	{
		//std::cout<<"column: "<<i<<"\n";
		if(GetNumber(CurrentRow, i)-1e-15 <= (T)0&&GetNumber(CurrentRow, i)+1e-15 >= (T)0)
		{
			//max=(T)0;
			//MaxRow=CurrentRow;
			for(j = CurrentRow+1; j < row; j++)
			{
				/*if((temp=abs(GetNumber(j, i))) > max)
				{
					max = temp;
					MaxRow=j;
				}*/
				if(GetNumber(j, i)-1e-15>=0||GetNumber(j, i)+1e-15<=0)
					break;
			}
			if(j == row)
				continue;
			else
			{
				RowChange(CurrentRow, j);
				Mirror.RowChange(CurrentRow, j);
				result = 0 - result;
				//std::cout<<"Rowchange: "<<CurrentRow<<" "<<j<<"\n";
			}
		}
		for(j = CurrentRow + 1; j < row; j++)
		{
			if(GetNumber(j, i) != (T)0)
			{
				/*std::cout<<"("<<j<<" "<<i<<")="<<GetNumber(j, i)<<"\n";
				std::cout<<"("<<CurrentRow<<" "<<i<<")="<<GetNumber(CurrentRow, i)<<"\n";
				std::cout<<"Add:"<<j<<" "<<CurrentRow<<" "<<(T)0-GetNumber(j, i)/GetNumber(CurrentRow, i)<<"\n";*/
				RowAdd(j, CurrentRow, temp=((T)0-GetNumber(j, i))/GetNumber(CurrentRow, i));
				Mirror.RowAdd(j, CurrentRow, temp);
			}
		}
		CurrentRow++;
	}
	//std::cout<<"Mirror:\n"<<Mirror;
	return result;
}

template<class T>
int MyMatrix<T>::GauessElimination()
{
	MyMatrix<T> A(row, row);
	return this->GauessElimination(A);
}

template<class T> class MySquareMatrix;

template<class T>
class MySquareMatrix: public MyMatrix<T>{
public:
	MySquareMatrix(int r, T* p): MyMatrix<T>(r, r, p) {}
	MySquareMatrix(int r): MyMatrix<T>(r, r){}
	MySquareMatrix<T>(const MySquareMatrix<T>& A): MyMatrix<T>(A) {}
	T determination() const;//���㷽�������ʽ 
	bool inverse(MySquareMatrix<T>& t) const;
	bool LU(MySquareMatrix<T>& L, MySquareMatrix<T>& U) const;
	bool sqrtLU(MySquareMatrix<T>& L, MySquareMatrix<T>& U) const;
	MyMatrix<T> SolveEqu(MyMatrix<T> b, int lu = 0) const;//��LU�ֽⷨ��� 
	MyMatrix<T> Jacobi(MyMatrix<T> b, int iter_number) const;//��Jacobi�������
	MyMatrix<T> Jacobi(MyMatrix<T> b, T err) const;
	MyMatrix<T> Gauss_Seidel(MyMatrix<T> b, int iter_number) const;//��GS�������
	MyMatrix<T> Gauss_Seidel(MyMatrix<T> b,T err) const;
	MyMatrix<T> SOR(MyMatrix<T> b, int iter_number, T omega) const;//��SOR�������
	MyMatrix<T> SOR(MyMatrix<T> b, T err, T omega) const;
	void tranverse();//ת��
	void rotate();//��ʱ����ת180�� 
	T infty_norm();
	T infty_cond();
	int GetRow() const
	{
		return (MyMatrix<T>::GetRow());
	} 
};

template<class T>
T MySquareMatrix<T>::infty_norm()
{
	T max = (T)0;
	int r = GetRow();
	for(int i = 0; i < r; i++)
	{
		T temp = (T)0;
		for(int j = 0; j < r; j++)
		{
			T element = MyMatrix<T>::GetNumber(i, j);
			if(element >= (T)0)
			{
				temp += element;
			}
			else
			{
				temp -= element;
			}
		}
		if(temp > max)
		{
			max = temp;
		}
	}
	return max;
}

template<class T>
T MySquareMatrix<T>::infty_cond()
{
	int r = GetRow();
	MySquareMatrix<T> B(r);
	if(!this->inverse(B))
	  return -1;
	//std::cout<<"here!"<<std::endl;
	return this->infty_norm() * B.infty_norm();
}

template<class T>
void MySquareMatrix<T>::tranverse()
{
	int i, j;
	T temp;
	int r = GetRow();
	for(i = 0; i < r; i++)
	{
		for( j = i + 1; j < r; j++)
		{
			temp = MyMatrix<T>::GetNumber(i, j);
			MyMatrix<T>::SetNumber(i, j, MyMatrix<T>::GetNumber(j, i));
			MyMatrix<T>::SetNumber(j, i, temp);
		}
	}
}

template<class T>
void MySquareMatrix<T>::rotate()
{
	int i, j;
	T temp;
	int r = GetRow();
	MySquareMatrix<T> A = *this;
	for(i = 0; i < r; i++)
	{
		for( j = 0; j < r; j++)
		{
			MyMatrix<T>::SetNumber(i, j, A.GetNumber(r - 1 - i, r - 1 - j));
		}
	}
}
template<class T>
T MySquareMatrix<T>::determination() const
{
	T id;
	MySquareMatrix<T> A = *this;
	id = (T)A.GauessElimination();
	int i;
	for(i = 0; i < A.GetRow(); i++)
	{
		id = id * A.GetNumber(i, i);
	}
	return id;
}

template<class T>
bool MySquareMatrix<T>::inverse(MySquareMatrix<T>& A) const
{
	int i,r;
	float temp;
	r=this->GetRow();
	for(i = 0; i < r; i++)
	{
		A.SetNumber(i, i, (T)1);
	}
	MySquareMatrix<T> B = *this;
	B.GauessElimination(A);
	temp = B.GetNumber(r-1,r-1);
	if(temp<=(T)1e-15&&temp>=-(T)1e-15)
		return false;
	//std::cout<<"B: \n"<<B;
	//std::cout<<"A: \n"<<A;
	B.rotate();
	A.rotate();
	B.GauessElimination(A);
	B.rotate();
	A.rotate();
	//std::cout<<"B: \n"<<B;
	//std::cout<<"A: \n"<<A;
	for(i = 0; i < r; i++)
	{
		A.RowMutiply(i, (T)1/B.GetNumber(i, i));
	}
	return true;
}

template<class T>
bool MySquareMatrix<T>::LU(MySquareMatrix<T>& L, MySquareMatrix<T>& U) const
{
	int row = this->GetRow();
	for(int k = 0; k < row; k++)
		U.SetNumber(0, k, this->GetNumber(0, k));
	for(int k = 0; k < row; k++)
	{
		L.SetNumber(k, 0, this->GetNumber(k, 0)/this->GetNumber(0, 0));
	}
	for(int r = 1; r < row; r++)
	{
		L.SetNumber(r, r, (T)1);
		for(int i = r; i < row; i++)
		{
			T sum = this->GetNumber(r, i);
			for(int k = 0; k < r; k++)
			{
				sum -= (L.GetNumber(r, k) * U.GetNumber(k, i));
			}
			U.SetNumber(r, i, sum);
		}
		for(int i = r+1; i < row; i++)
		{
			T sum = (T)0;
			for(int k = 0; k < r; k++)
			{
				sum += (L.GetNumber(i, k) * U.GetNumber(k, r));
			}
			if(U.GetNumber(r, r) == (T)0)
			{
				return false;
			}
			L.SetNumber(i, r, (this->GetNumber(i, r) - sum)/U.GetNumber(r, r));
		}
	}
	return true;
}

template<class T>
bool MySquareMatrix<T>::sqrtLU(MySquareMatrix<T>& L, MySquareMatrix<T>& U) const
{
	//������������ʹ��
	int row = this->GetRow();
	for(int j = 0; j < row; j++)
	{
		T temp_sum = (T)0;
		for(int k = 0; k < j; k++)
		{
			temp_sum += (pow(L.GetNumber(j, k), 2));
		}
		//cout<<this->GetNumber(j, j)<<" "<<temp_sum<<endl;
		T result1 = this->GetNumber(j, j) - temp_sum;
		if(result1 >= 0)
			L.SetNumber(j, j, sqrt(result1));
		else
			L.SetNumber(j, j, 0);
		for(int i = j + 1; i < row; i++)
		{
			T sum = (T)0;
			for(int k = 0; k < j; k++)
			{
				sum += L.GetNumber(i, k) * L.GetNumber(j, k);
			}
			L.SetNumber(i, j, (this->GetNumber(i, j) - sum)/L.GetNumber(j, j));
		}
	}
	U = L;
	U.tranverse();
	return true;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::Jacobi(MyMatrix<T> b, int iter_number) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	MyMatrix<T> next_x(row, 1);
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; t < iter_number; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			next_x.SetNumber(i, 0, temp);
		}
		current_x = next_x;
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::Jacobi(MyMatrix<T> b, T err) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	MyMatrix<T> next_x(row, 1);
	MyMatrix<T> difference(row, 1);
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; ; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			next_x.SetNumber(i, 0, temp);
		}
		difference = current_x - next_x;
		current_x = next_x;
		if(difference.infty_norm() < err)
		{
			cout<<"Jacobi iteration number: "<<t<<endl;
			break;
		}
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::Gauss_Seidel(MyMatrix<T> b, int iter_number) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; t < iter_number; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			current_x.SetNumber(i, 0, temp);
		}
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::Gauss_Seidel(MyMatrix<T> b, T err) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	T maximum_err = (T)0;
	T current_err = (T)0;
	int max_iter = 100;
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; t < max_iter; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			current_err = current_x.GetNumber(i, 0) - temp;
			current_err = current_err > (T)0 ? current_err : -current_err;
			current_x.SetNumber(i, 0, temp);
			if(current_err - maximum_err > (T)0)
			{
				maximum_err = current_err;
			}
		}
		if(maximum_err < err)
		{
			//cout << "GS iteration number: "<<t<<endl;
			break;
		}
		maximum_err = (T)0;
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::SOR(MyMatrix<T> b, int iter_number, T omega) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; t < iter_number; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			temp -= current_x.GetNumber(i, 0);
			temp *= omega;
			temp += current_x.GetNumber(i, 0);
			current_x.SetNumber(i, 0, temp);
		}
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::SOR(MyMatrix<T> b, T err, T omega) const
{
	int row = this->GetRow();
	MyMatrix<T> current_x(row, 1);
	T maximum_err = (T)0;
	T current_err = (T)0;
	int max_iter = 100000;
	for(int i = 0; i < row; i++)
	{
		current_x.SetNumber(i, 0, (T)0);
	}
	for(int t = 0; t < max_iter; t++)
	{
		for(int i = 0; i < row; i++)
		{
			T temp = b.GetNumber(i, 0);
			for(int j = 0; j < row; j++)
			{
				if(j != i)
					temp -= (this->GetNumber(i, j) * current_x.GetNumber(j, 0));
			}
			temp /= this->GetNumber(i, i);
			temp -= current_x.GetNumber(i, 0);
			temp *= omega;
			current_err = temp;
			temp += current_x.GetNumber(i, 0);
			current_err = current_err > 0 ? current_err : -current_err;
			current_x.SetNumber(i, 0, temp);
			if(current_err > maximum_err)
			{
				maximum_err = current_err;
			}
		}
		if(maximum_err < err)
		{
			cout<<"SOR iteration number: "<<t<<endl;
			break;
		}
		maximum_err = (T)0;
	}
	return current_x;
}

template<class T>
MyMatrix<T> MySquareMatrix<T>::SolveEqu(MyMatrix<T> b, int lu) const
{
	int row = this->GetRow();
	MySquareMatrix<T> L(row);
	MySquareMatrix<T> U(row);
	if(lu == 0)
		this->LU(L, U);
	else if(lu == 1)
	{
		this->sqrtLU(L, U);
	}
	vector<T> y;
	y.push_back(b.GetNumber(0, 0)/L.GetNumber(0, 0));
	for(int i = 1; i < row; i++)
	{
		//T temp_result = (T)0;
		y.push_back(b.GetNumber(i, 0));
		for(int j = 0; j < i; j++)
		{
			y[i] -= (L.GetNumber(i, j) * y[j]);
		}
		y[i] /= L.GetNumber(i, i);
	}
	/*for(int i = 0; i < row; i++)
	{
		std::cout<<y[i]<<" ";
	}*/
	T* x = new T[row];
	x[row - 1] = y[row - 1]/ U.GetNumber(row - 1, row - 1);
	for(int i = row - 2; i >= 0; i--)
	{
		x[i] = y[i];
		for(int j = i + 1; j < row; j++)
		{
			x[i] -= (U.GetNumber(i, j) * x[j]);
		}
		x[i] /= U.GetNumber(i, i);
		/*x[i] = (T)0;
		for(int j = i + 1; j < row; j++)
		{
			x[i] += (U.GetNumber(i, j) * x[j]);
		}
		x[i] = (y[i] - x[i])/U.GetNumber(i, i);*/
	}
	MyMatrix<T> result(row, 1);
	for(int i = 0; i < row; i++)
	{
		result.SetNumber(i, 0, x[i]);
	}
	delete[] x;
	return result;
}

const int DefaultSize = 20000;

template<class T>class SparseMatrix;
template<class T> std::ostream& operator<<(ostream& out,const SparseMatrix<T>& M);
template<class T> std::istream& operator>>(istream& in, SparseMatrix<T>& M);
template<class T>
struct Trituple{
	int row, col;
	T value;
	Trituple<T>& operator=(Trituple<T>& x)
	{
		row = x.row;
		col = x.col;
		value = x.value;
	}
};

template<class T>
class SparseMatrix{
	friend ostream& operator<< <>(ostream& out,const SparseMatrix<T>& M);
	friend istream& operator>> <>(istream& in, SparseMatrix<T>& M);
	//friend SparseMatrix<T>& operator=(SparseMatrix<T>& x);
	
	int Rows, Cols, Terms;
	Trituple<T>* smArray;
	int maxTerms;
	public:
		SparseMatrix(int r = 0, int c = 0, int t = 0,int maxSz = DefaultSize);
		SparseMatrix<T>(SparseMatrix<T>& x)//���ƹ��캯��
		{
			//cout<<"&\n";
			Rows = x.Rows; Cols = x.Cols; Terms = x.Terms; maxTerms =x.maxTerms;
			smArray = new Trituple<T>[maxTerms];
			for(int i = 0; i < Terms; i++)
			{
				this->smArray[i] = x.smArray[i];
			}
		}
		const SparseMatrix<T>& operator=(const SparseMatrix<T> &x)
		{
			Rows = x.Rows; Cols = x.Cols; Terms = x.Terms; maxTerms =x.maxTerms;
			delete[] smArray;
			smArray = new Trituple<T>[maxTerms];
			for(int i = 0; i < Terms; i++)
			{
				smArray[i] = x.smArray[i];
			}
			//std::cout<<"="<<endl;
			return *this;
		}
		/*SparseMatrix<T>(SparseMatrix<T>&& con)
		{
			Rows = con.Rows; Cols = con.Cols; Terms = con.Terms; maxTerms =con.maxTerms;
			smArray = con.smArray;
			con.smArray = NULL;
		}*/
		~SparseMatrix()
		{
			delete[] smArray;
		} 
		SparseMatrix<T> Transpose() const;
		SparseMatrix<T> Add(const SparseMatrix<T>& b) const;
		SparseMatrix<T> Multiply(const SparseMatrix<T>& b) const;
		SparseMatrix<T> Solve(const SparseMatrix<T>& b) const;//�ṩ����ĳ�ʼֵx��Ȼ���õ������ⷽ��
		void Set(int i, int r, int c, T v)//�� smArray[i]��Ϊr,c,v 
		{
			smArray[i].row = r;
			smArray[i].col = c;
			smArray[i].value = v;
		}
		T Get(int i)//ȡ��smArray[i]��value
		{
			return smArray[i].value;
		 }
		bool SolveAble()
		{
			int i;
			int t[Rows];
			for(i = 0; i < Rows; i++)
			{
				t[i] = 0;
			}
			for(i = 0; i < Terms; i++)
			{
				if(smArray[i].col == smArray[i].row)
				{
					t[smArray[i].col] = 1;
				}
			}
			for(i = 0; i < Rows; i++)
			{
				if(t[i] == 0)
					return false;
			}
			return true;
		 } 
};

template<class T>
SparseMatrix<T>::SparseMatrix(int r, int c, int t, int maxSz):Rows(r), Cols(c), Terms(t), maxTerms(maxSz)
{
	if(maxSz < 1)
	{
		cerr<<"�����ʼ������\n";
		exit(1);
	}
	if(Terms > maxTerms)
	{
		cerr<<"Terms overflow!\n";
	}
	smArray = new Trituple<T>[maxTerms];
	//Rows = Cols = Terms = 0;
}


template<class T>
ostream& operator<<(ostream& out,const SparseMatrix<T>& M)
{
	int i, j;
	int temp = 0;
	Trituple<T> Current = M.smArray[temp];
	for(i = 0; i < M.Rows; i++)
	{
		for(j = 0; j < M.Cols; j++)
		{
			if(Current.col != j||Current.row != i)
			{
				std::cout<<0<<' ';
			}
			else
			{
				std::cout<<Current.value<<' ';
				temp++;
				if(temp >= M.Terms)
				{
					Current.row = -1;
				}
				else
				{
					Current = M.smArray[temp];
				}
			}
		}
		cout<<endl;
	}
	cout<<"M.Terms:"<<M.Terms<<endl;
}

int comp(const void* c, const void* d)
{
	Trituple<double>* a = (Trituple<double>*)c;
	Trituple<double>* b = (Trituple<double>*)d;
	if(a->row < b->row)
	{
		return -1;
	}
	else if(a->row > b->row)
	{
		return 1;
	}
	else
	{
		if(a->col < b->col)
		{
			return -1;
		}
		else if(a->col > b->col)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
}
//����Ҫ�����ʹ��cin,�������ù��캯���Թ�ģ�������������г�ʼ��
//��������У������С��С���ֵ�����У�������Terms���� 
template<class T>
istream& operator>>(istream& in, SparseMatrix<T>& M)
{
	for(int i = 0; i < M.Terms; i++)
	{
		//in>>M.smArray[i].row>>M.smArray[i].col>>M.smArray[i].value;
		int t, j;
		T v;
		in>>t>>j>>v;
		M.smArray[i].row = t-1;
		M.smArray[i].col = j-1;
		M.smArray[i].value = v;
	}
	qsort(M.smArray, M.Terms, sizeof(Trituple<double>), comp);
	return in;
}


template<class T>
SparseMatrix<T> SparseMatrix<T>::Transpose() const
{
	int* rowSize = new int[Cols];//ͳ��ÿ�з���Ԫ�صĸ���
	int* rowStart = new int[Cols];//Ԥ��ת�ú��������λ��
	SparseMatrix<T> b(Cols, Rows, Terms, maxTerms);//���ת�ý��
	if(Terms > 0)
	{
		int i, j;
		for(i = 0; i < Cols; i++)rowSize[i] = 0;
		for(i = 0; i < Terms; i++)rowSize[smArray[i].col]++;
		rowStart[0] = 0;
		for(i = 1; i < Cols; i++)
			rowStart[i] = rowStart[i-1] + rowSize[i-1];
		for(i = 0; i < Terms; i++)
		{
			j = rowStart[smArray[i].col];
			b.smArray[j].row = smArray[i].col;
			b.smArray[j].col = smArray[i].row;
			b.smArray[j].value = smArray[i].value;
			rowStart[smArray[i].col]++;
		}
	}
	delete[] rowSize;
	delete[] rowStart;
	//std::cout<<b;
	return b;
}

template<class T>
SparseMatrix<T> SparseMatrix<T>::Add(const SparseMatrix<T>& b) const
{
	SparseMatrix<T> result(Rows, Cols);
	if(Rows != b.Rows || Cols != b.Cols)
	{
		cout<<"���һ�£��޷���ӣ�\n";
		return result;
	}
	int i = 0, j = 0, t = 0, index_a, index_b;
	while(i < Terms && j < b.Terms)
	{
		if(comp(&smArray[i], &b.smArray[j]) == -1)
		{
			result.smArray[t] = smArray[i];
			i++;
			t++;
		}
		else if(comp(&smArray[i], &b.smArray[j]) == 1)
		{
			result.smArray[t] = smArray[j];
			j++;
			t++;
		}
		else
		{
			result.smArray[t] = smArray[i];
			result.smArray[t].value = smArray[i].value + b.smArray[j].value;
			i++;
			j++;
			t++;
		}
	}
	for(; i < Terms; i++)
	{
		result.smArray[t] = smArray[i];
		t++;
	}
	for(; j < b.Terms; j++)
	{
		result.smArray[t] = b.smArray[j];
		t++;
	}
	result.Terms = t;
	return result;
}

template<class T>
SparseMatrix<T> SparseMatrix<T>::Multiply(const SparseMatrix<T>& b) const
{
	int i;
	SparseMatrix<T> result(Rows, b.Cols);
	SparseMatrix<T> c(b.Cols, b.Rows, b.Terms);
	if(Cols != b.Rows)
	{
		cerr<<"Incompatible matrices\n";
		return result;
	}
	if(Terms == maxTerms || b.Terms == maxTerms)
	{
		cerr<<"One additional space in a or b needed"<<endl;
		return result;
	}
	int* rowSize = new int[b.Rows];//����B���з���Ԫ�صĸ��� 
	int* rowStart = new int[b.Rows + 1];//����B��������Ԫ���п�ʼ��λ�� 
	T* temp = new T[b.Cols];
	int Current, lastInResult, RowA, ColA, ColB;
	for(i = 0; i < b.Rows; i++) 
		rowSize[i] = 0;
	for(i = 0; i < b.Terms; i++) 
		rowSize[b.smArray[i].row]++;
	rowStart[0] = 0;
	for(i = 1; i <= b.Rows; i++)
	{
		rowStart[i] = rowStart[i-1] + rowSize[i-1];
	}
	Current = 0;lastInResult = -1; 
	while(Current < Terms)
	{
		RowA = smArray[Current].row;
		for(i = 0; i < b.Cols; i++)temp[i] = 0;
		while(Current < Terms && smArray[Current].row == RowA)
		{
			ColA = smArray[Current].col;
			for(i = rowStart[ColA]; i < rowStart[ColA+1]; i++)
			{
				ColB = b.smArray[i].col;
				temp[ColB] += smArray[Current].value*b.smArray[i].value;
			}
			Current++;
		}
		for(i = 0; i < b.Cols; i++)
		{
			if(temp[i] - 0.001>= 0||temp[i]+0.001<=0)
			{
				lastInResult++;
				result.smArray[lastInResult].row = RowA;
				result.smArray[lastInResult].col = i;
				result.smArray[lastInResult].value = temp[i];
			}
		}
	}
	//std::cout<<"lastInResult:"<<lastInResult<<endl;
	result.Rows = Rows;result.Cols = b.Cols;
	result.Terms = lastInResult+1;
	//cout<<result.smArray[0].row<<" "<<result.smArray[0].value<<endl;
	delete[] rowSize;
	delete[] rowStart;
	delete[] temp;
	return result;
}

template<class T>
SparseMatrix<T> SparseMatrix<T>::Solve(const SparseMatrix<T>& b) const
{
	SparseMatrix<T> S(Rows, Cols, Rows);//��ȡ�ĶԽ������ 
	SparseMatrix<T> R(Rows, Cols);//��ȥ�Խ����ľ���ĸ�ֵ
	int i, j, Rstart, Sstart; 
	for(i = 0; i < Rows; i++)
	{
		S.smArray[i].row = i;
		S.smArray[i].col = i;
		S.smArray[i].value = (T)0.001;
	}
	//cout<<*this;
	Sstart = 0; Rstart = 0;
	for(i = 0; i < Terms; i++)
	{
		if(smArray[i].row == smArray[i].col)
		{
			/*S.smArray[Sstart].row = smArray[i].row;
			S.smArray[Sstart].col = smArray[i].col;*/
			S.smArray[smArray[i].row].value = smArray[i].value;
			Sstart++;
		}
		else
		{
			R.smArray[Rstart].row = smArray[i].row;
			R.smArray[Rstart].col = smArray[i].col;
			R.smArray[Rstart].value = -smArray[i].value;
			Rstart++;
		}
	}
	//cout<<"Sstart:"<<Sstart<<endl;
	R.Terms = Rstart+1;
	//std::cout<<"S\n"<<S;
	//std::cout<<"R\n"<<R;
	SparseMatrix<T> x(Rows, 1, Rows);
	SparseMatrix<T> temp(Rows, 1, Rows);
	for(i = 0; i < Rows; i++)
	{
		temp.smArray[i].row = i;
		temp.smArray[i].col = 0;
	}
	SparseMatrix<T> temp2(Rows, 1, Rows);
	//SparseMatrix<T> temp3(Rows, 1, Rows);
	for(i = 0; i < x.Terms; i++)
	{
		x.smArray[i].row = i;
		x.smArray[i].col = 0;
		x.smArray[i].value = b.smArray[i].value/S.smArray[i].value;
	}
	//cout<<"x\n"<<x;
	bool flag;
	T diff;
	for(i = 0; i < 500; i++)
	{
		temp2 = b.Add(R.Multiply(x));
		for(j = 0; j < Rows; j++)
		{
			temp.smArray[j].value = temp2.smArray[j].value/S.smArray[j].value;
			//cout<<temp.smArray[j].value<<" "<<temp.smArray[j].row<<" "<<temp.smArray[j].col<<endl;
		}
		//cout<<"temp.Terms"<<temp.Terms<<endl;
		//cout<<"temp\n"<<temp;
		flag = true;
		for(j = 0; j < Rows; j++)
		{
			diff = x.smArray[j].value-temp.smArray[j].value;
			if(diff-0.000001 > 0||diff+0.000001 < 0)
			{
				flag = false;
				break;
			}
		}
		if(flag == true)
		{
			break;
		}
		else
		{
			x = temp;
			//cout<<"x\n"<<x;
		}
	}
	return x;
}
