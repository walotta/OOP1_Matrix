#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>

using std::size_t;

namespace sjtu
{
	template <class T>
	class Matrix
	{

        template <class V, class U>
        friend auto operator*(const Matrix<V> &mat, const U &x);

        template <class V, class U>
        friend auto operator*(const U &x, const Matrix<V> &mat);

        template <class U, class V>
        friend auto operator*(const Matrix<U> &a, const Matrix<V> &b);

        template <class U, class V>
        friend auto operator+(const Matrix<U> &a, const Matrix<V> &b);

        template <class U, class V>
        friend auto operator-(const Matrix<U> &a, const Matrix<V> &b);

	private:
		// your private member variables here.
		size_t matrix_n=0;//行
		size_t matrix_m=0;//列
		size_t total_number=0;
		T* head= nullptr;//头指针
		
	public:
		Matrix() = default;
		
		Matrix(size_t n, size_t m, T _init = T())
		{
		    if(n<0||m<0)throw std::invalid_argument("CAN'T BUILD");
		    matrix_n=n;
		    matrix_m=m;
		    total_number=m*n;
            if(total_number==0)return;
            head=new T[total_number];
			for(size_t i=0;i<total_number;i++)
            {
                head[i]=_init;
            }
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
		    size_t n=sz.first;
		    size_t m=sz.second;

            if(n<0||m<0)throw std::invalid_argument("CAN'T BUILD");
            matrix_n=n;
            matrix_m=m;
            total_number=m*n;
            if(total_number==0)return;
            head=new T[total_number];
            for(size_t i=0;i<total_number;i++)
            {
                head[i]=_init;
            }
		}
		
		Matrix(const Matrix &o)
		{
			matrix_n=o.matrix_n;
			matrix_m=o.matrix_m;
			total_number=o.total_number;

            if(o.head== nullptr)return;
            head=new T[total_number];
			for(size_t i=0;i<total_number;i++)
            {
                head[i]=o.head[i];
            }
		}
		
		template <class U>
		Matrix(const Matrix<U> &o)
		{
            matrix_n=o.rowLength();
            matrix_m=o.columnLength();
            total_number=matrix_n*matrix_m;

            if(o.rowLength()== 0)return;
            head=new T[total_number];
            for(size_t i=0;i<total_number;i++)
            {
                head[i]=o(i/o.columnLength(),i%o.columnLength());
            }
		}
		
		Matrix &operator=(const Matrix &o)
		{
            clear();
            matrix_m=o.matrix_m;
            matrix_n=o.matrix_n;
            total_number=o.total_number;
            if(head!= nullptr)delete []head;
            if(o.total_number==0)return *this;
            head=new T[total_number];
            for(size_t i=0;i<total_number;i++)
            {
                head[i]=o.head[i];
            }
            return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
            clear();
            matrix_n=o.rowLength();
            matrix_m=o.columnLength();
            total_number=matrix_n*matrix_m;
            if(head!= nullptr)delete []head;
            if(total_number==0)return *this;
            head=new T[total_number];
            for(size_t i=0;i<total_number;i++)
            {
                head[i]=o(i/o.columnLength(),i%o.columnLength());
            }
            return *this;
		}
		
		Matrix(Matrix &&o) noexcept
		{
            matrix_m=o.matrix_m;
            matrix_n=o.matrix_n;
            total_number=o.total_number;
            head=o.head;
            o.head= nullptr;
		}
		
		Matrix &operator=(Matrix &&o) noexcept
		{
			clear();
            matrix_n=o.matrix_n;
            matrix_m=o.matrix_m;
            total_number=o.total_number;
            if(head!= nullptr)delete []head;
            head=o.head;
            o.head= nullptr;
            return *this;
		}
		
		~Matrix() {
		    matrix_n=0;
		    matrix_m=0;
		    total_number=0;
		    if(head!= nullptr)delete []head;
		}
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
		    matrix_n=il.size();
            auto it_row=il.begin();
            matrix_m=it_row->size();
            if(matrix_n<=0||matrix_m<=0)throw std::invalid_argument("BUILD FAILED");
            total_number=matrix_m*matrix_n;
            if(total_number==0)return;
            head=new T[total_number];
			for(size_t i=0;i<matrix_n;i++)
            {
                auto it=it_row->begin();
                for(size_t j=i*matrix_m;j<(i+1)*matrix_m;j++)
                {
                    head[j]=*it;
                    it++;
                }
                it_row++;
            }
		}
		
	public:
		size_t rowLength() const {
		    return matrix_n;
		}
		
		size_t columnLength() const {
            return matrix_m;
		}
		
		void resize(size_t _n, size_t _m, T _init = T())
		{
		    if(_n<=0||_m<=0)throw std::invalid_argument("INVALID RESIZE");
			if(_n*_m==total_number)
            {
			    matrix_m=_m;
			    matrix_n=_n;
                return;
            }else
            {
			    T* tem;
			    size_t new_total=_n*_m;
			    tem=new T[new_total];
                for(size_t i=0;i<new_total;i++)
                {
                    if(i<total_number)tem[i]=head[i];
                    else tem[i]=_init;
                }
                matrix_m=_m;
                matrix_n=_n;
                total_number=new_total;
                delete []head;
                head=tem;
                tem= nullptr;
                return;
            }
		}
		
		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
			size_t _n=sz.first;
			size_t _m=sz.second;

            if(_n<=0||_m<=0)throw std::invalid_argument("INVALID RESIZE");
            if(_n*_m==total_number)
            {
                matrix_m=_m;
                matrix_n=_n;
                return;
            }else
            {
                T* tem;
                size_t new_total=_n*_m;
                tem=new T[new_total];
                for(size_t i=0;i<new_total;i++)
                {
                    if(i<total_number)tem[i]=head[i];
                    else tem[i]=_init;
                }
                matrix_m=_m;
                matrix_n=_n;
                total_number=new_total;
                delete []head;
                head=tem;
                tem= nullptr;
                return;
            }
		}
		
		std::pair<size_t, size_t> size() const
		{
			std::pair<size_t,size_t> tem(matrix_n,matrix_m);
			return tem;
		};
		
		void clear()
		{
			matrix_m=0;
			matrix_n=0;
			total_number=0;
			if(head!= nullptr)delete []head;
			head= nullptr;
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
		    //todo why
		    if(i>=matrix_n||j>=matrix_m||i<0||j<0)throw std::invalid_argument("TOO LARGE POS");
            return head[i*matrix_m+j];
		}
		
		T &operator()(size_t i, size_t j)
		{
			//todo why
            if(i>=matrix_n||j>=matrix_m||i<0||j<0)throw std::invalid_argument("INVALID POS");
            return head[i*matrix_m+j];
		}
		
		Matrix<T> row(size_t i) const
		{
		    if(i>=matrix_n||i<0)throw std::invalid_argument("INVALID ROW NUMBER");
			Matrix tem(1,matrix_m);
			for(int j=0;j<matrix_m;j++)
            {
                tem.head[j]=head[j+i*matrix_m];
            }
            return tem;
		}
		
		Matrix<T> column(size_t i) const
		{
			if(i>=matrix_m||i<0)throw std::invalid_argument("INVALID COLUMN NUMBER");
			Matrix tem(matrix_n,1);
			for(size_t j=0;j<matrix_n;j++)
            {
			    tem.head[j]=head[j*matrix_m+i];
            }
            return tem;
		}
		
		
	public:
		template <class U>
		bool operator==(const Matrix<U> &o) const
		{
            if(size()!=o.size())return false;
            for(size_t i=0;i<total_number;i++)
            {
                if(head[i]!=o(i/matrix_m,i%matrix_m))return false;
            }
            return true;
		}
		
		template <class U>
		bool operator!=(const Matrix<U> &o) const
		{
            return !(*this==o);
		}
		
		Matrix operator-() const
		{
		    for(size_t i=0;i<total_number;i++)
            {
		        head[i]=-head[i];
            }
            return *this;
		}
		
		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
			if(size()!=o.size())throw std::invalid_argument("INVALID ADD");
            *this=*this+o;
            return *this;
		}
		
		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
            if(size()!=o.size())throw std::invalid_argument("INVALID SUB");
            *this=*this-o;
            return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
            *this=(*this)*x;
            return *this;
		}
		
		Matrix tran() const
		{
			Matrix TEM(matrix_m,matrix_n);
			for(size_t i=0;i<total_number;i++)
            {
			    TEM(i/matrix_n,i%matrix_n)=(*this)(i%matrix_n,i/matrix_n);
            }
			return TEM;
		}
		
	public: // iterator
	    //todo iterator
		class iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			iterator() = default;
			
			iterator(const iterator &) = default;
			
			iterator &operator=(const iterator &) = default;

			~iterator()
            {
			    n_it=0;
			    m_it=0;
                if(private_matrix&&now_matrix!= nullptr)delete []now_matrix;
            }

			friend iterator Matrix::begin();
			friend iterator Matrix::end();
			friend std::pair<iterator, iterator> Matrix::subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r);
			
		private:
            size_type n_it;
            size_type m_it;
            Matrix* now_matrix;
            bool private_matrix=false;
			
		public:
            size_type id() const
            {
                return m_it+n_it*now_matrix->columnLength();
            }

			difference_type operator-(const iterator &o)
			{
                if(now_matrix!=o.now_matrix)throw std::invalid_argument("NOT SAME MATRIX");
				return this->id()-o.id();
			}
			
			iterator &operator+=(difference_type offset)
			{
                if(n_it*now_matrix->columnLength()+m_it+offset>now_matrix->columnLength()*now_matrix->rowLength())throw std::invalid_argument("TOO LARGE OFFSET");
				n_it+=(m_it+offset)/now_matrix->columnLength();
				m_it=(m_it+offset)%now_matrix->columnLength();
                return *this;
			}
			
			iterator operator+(difference_type offset) const
			{
                if(n_it*now_matrix->columnLength()+m_it+offset>now_matrix->columnLength()*now_matrix->rowLength())throw std::invalid_argument("TOO LARGE OFFSET");
                iterator tem_it;
                tem_it.now_matrix=now_matrix;
                tem_it.n_it=(m_it+offset)/now_matrix->columnLength();
                tem_it.m_it=(m_it+offset)%now_matrix->columnLength();
                return tem_it;
			}
			
			iterator &operator-=(difference_type offset)
			{
                if(n_it*now_matrix->columnLength()+m_it-offset<0)throw std::invalid_argument("TOO LARGE OFFSET");
                n_it=(n_it*now_matrix->columnLength()+m_it-offset)/now_matrix->columnLength();
                m_it=(n_it*now_matrix->columnLength()+m_it-offset)%now_matrix->columnLength();
                return *this;
			}
			
			iterator operator-(difference_type offset) const
			{
                if(n_it*now_matrix->columnLength()+m_it-offset<0)throw std::invalid_argument("TOO LARGE OFFSET");
                iterator tem_it;
                tem_it.now_matrix=now_matrix;
                tem_it.n_it=(n_it*now_matrix->columnLength()+m_it-offset)/now_matrix->columnLength();
                tem_it.m_it=(n_it*now_matrix->columnLength()+m_it-offset)%now_matrix->columnLength();
                return tem_it;
			}
			
			iterator &operator++()
			{
				*this=*this+1;
				return *this;
			}
			
			iterator operator++(int)
			{
				iterator tem_it;
				tem_it=*this;
				*this=*this+1;
                return tem_it;
			}
			
			iterator &operator--()
			{
                *this=*this-1;
                return *this;
			}
			
			iterator operator--(int)
			{
                iterator tem_it;
                tem_it=*this;
                *this=*this-1;
                return tem_it;
			}
			
			reference operator*() const
			{
			    return (*now_matrix)(n_it,m_it);
			}
			
			pointer operator->() const
			{
                return &((*now_matrix)(n_it,m_it));
			}
			
			bool operator==(const iterator &o) const
			{
                if(n_it!=o.n_it||m_it!=o.m_it||now_matrix!=o.now_matrix)return false;
                else return true;
			}
			
			bool operator!=(const iterator &o) const
			{
				return !(*this==o);
			}
		};
		
		iterator begin()
		{
			iterator tem_it;
            tem_it.now_matrix=this;
            tem_it.m_it=0;
            tem_it.n_it=0;
            return tem_it;
		}
		
		iterator end()
		{
            iterator tem_it;
            tem_it.now_matrix=this;
            tem_it.m_it=matrix_m-1;
            tem_it.n_it=matrix_n-1;
            return tem_it;
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r)
		{
			//todo may wrong
			size_t n_length;
			size_t m_length;
			n_length=r.first-l.first;
			m_length=r.second-l.second;
			if(n_length<=0||m_length<=0||n_length>matrix_n||m_length>matrix_m)throw std::invalid_argument("INVALID POS");
            std::pair<iterator, iterator> pair_return;
			iterator it_begin;
			iterator it_end;
			it_begin.m_it=0;
			it_begin.n_it=0;
			it_begin.private_matrix=true;
			it_end.n_it=n_length-1;
			it_end.m_it=m_length-1;
			it_end.private_matrix=true;
			it_begin.now_matrix=new Matrix(n_length,m_length);
			it_end.now_matrix=new Matrix(n_length,m_length);
            for(size_t i=0;i<it_begin.now_matrix->total_number;i++)
            {
                it_begin.now_matrix->head[i]=(*this)(l.first+i/m_length,l.second+i%m_length);
                it_end.now_matrix->head[i]=(*this)(l.first+i/m_length,l.second+i%m_length);
            }
            pair_return.first=it_begin;
            pair_return.second=it_end;
			return pair_return;
        }
	};

}

//
namespace sjtu
{
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)
	{
	    T t_ex;
	    U u_ex;
	    Matrix<decltype(t_ex*u_ex)> return_matrix(mat);
		for(size_t i=0;i<mat.columnLength()*mat.rowLength();i++)
        {
		    return_matrix.head[i]=x*return_matrix.head[i];
        }
		return return_matrix;
	}
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)
	{
        T t_ex;
        U u_ex;
        Matrix<decltype(t_ex*u_ex)> return_matrix(mat);
        for(size_t i=0;i<mat.columnLength()*mat.rowLength();i++)
        {
            return_matrix.head[i]=x*return_matrix.head[i];
        }
        return return_matrix;
	}
	
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)
	{
        if(a.columnLength()!=b.rowLength())throw std::invalid_argument("MATRIX NOT FIT");
        U u_ex;
        V v_ex;
        Matrix<decltype(u_ex*v_ex)> return_matrix(a.matrix_n,b.matrix_m);
        for(size_t y=0;y<return_matrix.matrix_n;y++)
        {
            for(size_t x=0;x<return_matrix.matrix_m;x++)
            {
                size_t tem_sum=0;
                for(size_t i=0;i<a.matrix_m;i++)
                {
                    tem_sum+=a(x,i)*b(i,y);
                }
                return_matrix(x,y)=tem_sum;
            }
        }
        return return_matrix;
	}
	
	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)
	{
		if(a.columnLength()!=b.columnLength()||a.rowLength()!=b.rowLength())throw std::invalid_argument("MATRIX NOT FIT");
		U tem_a;
		V tem_b;
        Matrix<decltype(tem_a+tem_b)> return_matrix(a);
        for(size_t i=0;i<a.total_number;i++)
        {
            return_matrix.head[i]=a.head[i]+b.head[i];
        }
        return return_matrix;
	}
	
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)
	{
        if(a.columnLength()!=b.columnLength()||a.rowLength()!=b.rowLength())throw std::invalid_argument("MATRIX NOT FIT");
        U tem_a;
        V tem_b;
        Matrix<decltype(tem_a+tem_b)> return_matrix(a);
        for(size_t i=0;i<a.total_number;i++)
        {
            return_matrix.head[i]=a.head[i]-b.head[i];
        }
        return return_matrix;
	}
	
}

#endif //SJTU_MATRIX_HPP

