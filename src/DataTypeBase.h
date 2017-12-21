#ifndef GFC_DATATYPEBASE_HPP
#define GFC_DATATYPEBASE_HPP

#include "GException.h"

namespace gfc
{
	
	//所有的数据类型的基类
	class DataTypeBase
	{
		
	public:
		/// virtual desctuctor
		
		virtual ~DataTypeBase(void) {};
		
		virtual void dump(std::ostream& s) const {}
		
		//派生类必须实现该函数getClassName()
		virtual GString getClassName(void) const = 0;
				
	private:
		GString m_typeName;   //该数据类型的名称
				
	};
	
	
}


#endif