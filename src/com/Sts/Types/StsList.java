
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.Reflect.*;
import com.Sts.Utilities.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import java.lang.reflect.*;
import java.util.*;

/** A general list class implemented by StsObjectList. Designed to be an efficient
 *  and growable index-accessible list.  Methods were named in keeping with
 *  Collection class style in JDK 1.2.  This can be easily substituted with
 *  a standard Java class when efficient/simple ones are available with 1.2
 */
public class StsList implements StsSerializable
{
    public Object[] list = null;
    protected int size = 0;
    protected int incSize;
    protected int currentPosition = 0;

    static public final int NO_MATCH = StsParameters.NO_MATCH;

    public StsList()
    {
    }

    public StsList(int initSize, int incSize)
    {
        this.incSize = incSize;
        initList(initSize);
    }

   	public StsList(int initSize)
    {
        this(initSize, 1);
    }

    public StsList(StsList list)
    {
        this.list = getList();
        this.size = list.size;
        this.incSize = list.incSize;
        this.currentPosition = list.currentPosition;
    }

      private void initList(int initSize)
    {
    	list = new Object[initSize];
    }

    public Object[] getList() { return list; }

    public Object[] getTrimmedList()
    {
        if(size == list.length) return list;
        Object[] trimmedList = new Object[size];
        System.arraycopy(list, 0, trimmedList, 0, size);
        return trimmedList;
    }

	public Object getCastList(Class c)
	{
		return StsMath.arraycopy(list, size, c);
	}

    public void setList(Object[] list) { this.list = list; }

    public Object copyArrayList()
    {
        return StsMath.arraycopy(list, size);
    }

    public Object copyArrayList(Class componentType)
     {
         return StsMath.arraycopy(list, size, componentType);
     }

    public Object getElement(int index)
    {

        if(index >= size || index < 0) return null; // object deleted or never existed
        if(!checkIndex(index, "StsList.getElement(int)")) return null;
        return list[index];
    }

    public void setElement(Object obj, int index)
    {
        if(!checkIndex(index, "StsList.getElement(int)")) return;
        list[index] = obj;
        if(size < index+1) size = index + 1;
    }

    private boolean checkIndex(int index, String methodDescription)
    {
    	if(index >= list.length)
        {
            StsException.outputException( new StsException(StsException.WARNING,
              "StsList.checkIndex() failed.", "index: " + index + " > list.length: " + list.length));
            return false;
        }
        else if(index < 0)
        {
            StsException.outputException(  new StsException(StsException.WARNING,
              "StsList.checkIndex() failed.", "index: " + index + " < 0"));
            return false;
        }
        return true;
    }

    public void setSize(int size) { this.size = size; }

    public int add(Object obj)
    {
    	int index;

        index = size;
        if (list.length <= size) resize(list.length);
        size++;

        try { list[index] = obj; }
        catch(Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return NO_MATCH;
        }
        return index;
    }

    public void add(Object[] objs)
    {
        if (objs==null) return;
        for (int i=0; i<objs.length; i++) add(objs[i]);
    }

    public boolean delete(int index)
    {
      	if(!checkIndex(index, "StsList.delete(int)")) return false;
        if(index >= size) return false; // object already deleted.

        size--;

        for(int n = index; n < size; n++)
          list[n] = list[n+1];

        list[size] = null;
        return true;
    }

    public void deleteLast()
    {
        if(size == 0) return;
        delete(size-1);
    }

    public int delete(Object obj)
    {
        for( int i=0; i<list.length; i++ )
        {
            if( list[i] == obj )
            {
                delete(i);
                return i;
            }
        }
        return NO_MATCH;
    }

 	public void deleteList(StsList deleteList)
    {
        int deleteSize = deleteList.getSize();

        for(int n = 0; n < deleteSize; n++)
         	delete((Object)deleteList.list[n]);
    }

    public void clear()
    {
        for(int i = 0; i < size; i++)
            list[i] = null;
        size = 0;
    }

 	public void addList(StsList extraList)
    {
        if(extraList == null) return;

        int extraSize = extraList.getSize();
      	if (list.length < size + extraSize)
            resize(size + extraSize);

        for(int n = 0; n < extraSize; n++)
         	add((Object)extraList.list[n]);
    }

 	public void addListReverse(StsList extraList)
    {
        if(extraList == null) return;

        int extraSize = extraList.getSize();
      	if (list.length < size + extraSize)
            resize(size + extraSize);

        for(int n = extraSize-1; n >= 0; n--)
         	add((Object)extraList.list[n]);
    }

    public void resize(int newSize)
    {
        Object[] oldList = list;
        list = new Object[newSize + incSize];
        System.arraycopy(oldList, 0, list, 0, oldList.length);
   	}

	public int insertBefore(int index, Object newObj)
    {
        if (list.length <= size) resize(list.length);

        // don't resize if first element beyond end of list is null
        if(list[size-1] != null) size++;

        for(int n = size-1; n > index; n--)
        	list[n] = list[n-1];

        list[index] = newObj;
        return index;
    }

	public int insertBefore(Object obj, Object newObj)
    {
        int index = getIndex(obj);

        if(index == NO_MATCH)
        {
            StsException.outputException( new StsException(StsException.WARNING,
                "StsList.insertBefore(Object, newObject) failed.",
                "Couldn't find matching Object in list."));
        }
        return insertBefore(index, newObj);
    }

    public int insertAfter(int index, Object newObj)
    {
        if (list.length <= size) resize(list.length);
        size++;

        for(int n = list.length-1; n > index+1; n--)
        	list[n] = list[n-1];

        list[index+1] = newObj;
        return index+1;
    }

    public void insertAfter(int index, StsList extraList)
    {
        if(extraList == null) return;
        int extraSize = extraList.getSize();
        if(extraSize == 0) return;

        Object[] extraListArray = extraList.getList();

        int newSize = size + extraSize;
      	if (list.length < newSize) resize(newSize);

        int nItemsToMove = size - index - 1;
        if(nItemsToMove > 0)
            System.arraycopy(list, index+1, list, index + extraSize + 1, nItemsToMove);
        System.arraycopy(extraListArray, 0, list, index+1, extraSize);

        setSize(newSize);
    }

    public void unInsert(Object obj)
    {
        int index = getIndex(obj);
        for(int i=index; i<list.length-1; i++)
        	list[i] = list[i+1];
        size--;
    }

	public int insertAfter(Object obj, Object newObj)
    {
        int index = getIndex(obj);

        if(index == NO_MATCH)
        {
            StsException.outputException( new StsException(StsException.WARNING,
                "StsList.StsList.insertAfter(Object, newObject) failed.",
                "Couldn't find matching Object in list."));
        }
        return insertAfter(index, newObj);
    }

    public int insertAt(int index, Object newObj)
    {
        if (list.length <= index+1)
        {
            resize(list.length);
            size = index+1;
        }
        list[index] = newObj;
        return index;
    }

    public int getIndex(Object object)
    {
        if (object==null) return NO_MATCH;
    	for (int n = 0; n < list.length; n++)
        	if (object.equals(list[n])) return n;

        return NO_MATCH;
    }


   	public int getSize() { return size; }
    public int getIncSize() { return incSize; }

	static public int getSize(StsList list)
	{
		if(list == null) return 0;
		else return list.getSize();
	}

    public Object getLast()
    {
/*
        if(size <= 0) return null;
        else return list[size-1];
*/
// should really use above...not sure why we are doing this. TJL 6/30/05
    	Object last = null;
        for( int i=list.length-1; i>=0; i-- )
        	if( list[i] != null )
            {
            	currentPosition = i;
                last = list[i];
                break;
            }
    	return last;
    }

    public Object getSecondToLast()
    {
        if(size <= 1) return null;
        return list[size-2];
    }

    public Object getFirst()
    {
/*
        if(size <= 0) return null;
        return list[0];
*/
    	Object first = null;
        for( int i=0; i<list.length; i++ )
        	if( list[i] != null )
            {
            	currentPosition = i;
                first = list[i];
                break;
            }

    	return first;
    }

    public Object getSecond()
    {
        if(size <= 1) return null;
        return list[1];
    }

    public Object getNext()
    {
    	Object next = null;
        currentPosition++;
		for( int i=currentPosition; i<list.length; i++ )
		if( list[i] != null )
        {
			currentPosition = i;
            next = list[i];
            break;
        }
        return next;
    }

    public void deleteAll()
    {
		if(list != null)
		{
			for (int n = 0; n < list.length; n++)
				list[n] = null;
		}
        size = 0;
    }

    public boolean contains(Object object)
    {
        if( object == null ) return false;

    	for( int i=0; i<list.length; i++ )
        {
        	if( list[i] != null )
            {
            	if( object.equals(list[i]) )
                	return true;
            }
        }
        return false;
    }

    public Object setCurrentPosition(int index)
    {
        if(!checkIndex(index, "StsList.setCurrentPosition(int)")) return null;
        currentPosition = index;
        return list[index];
    }
/*
    public void reverseOrder()
    {
        int n0 = 0;
        int n1 = size-1;

        while(n0 < n1)
        {
            Object item = list[n0];
            list[n0] = list[n1];
            list[n1] = item;
            n0++; n1--;
        }
    }
*/
    public void reverseOrder()
    {
        int n0, n1;

        if(list.length == 0) return;

        Object[] newList = new Object[list.length];

        for(n0 = 0, n1 = size-1; n0 < size; n0++, n1--)
            newList[n0] = list[n1];

        list = newList;
    }

 	public void forEach(String methodName, Object[] args)
    {
        if(list == null) return;

        int noElements = size;

        for(int n = 0; n < noElements; n++)
        {
            Object object = list[n];
		    if(object == null) continue;
            try
            {
                StsMethod method = new StsMethod(object, methodName, args);
                Method m = method.getMethod();
        		m.invoke(object, args);
            }
            catch(Exception e)
            {
        	    StsException.outputException("Method: " + methodName + " failed " +
                                             "for object of class: " +
                                             object.getClass().getName(),
                                             e, StsException.WARNING);
            }
        }
    }

 	public void forEach(String methodName)
    {
    	forEach(methodName, null);
    }

    /** If a single argument, make an array of length 1 for compatibility */
 	public void forEach(String methodName, Object arg)
    {
    	forEach(methodName, new Object[]{arg});
    }

	public void forEachInvoke(Method method, Object[] args)
   {
	   if(list == null) return;

	   int noElements = size;

	   for(int n = 0; n < noElements; n++)
	   {
		   Object object = list[n];
		   if(object == null) continue;
		   try
		   {
			   method.invoke(object, args);
		   }
		   catch(Exception e)
		   {
			   StsException.outputException("Method: " + method.getName() + " failed " +
											"for object of class: " +
											object.getClass().getName(),
											e, StsException.WARNING);
		   }
	   }
   }

    public void setValue(String valueName, Class valueClass, Object value)
    {
        if(list == null || size == 0) return;

        Class objectClass = null;
        Method setMethod;
        try
        {
            objectClass = list[0].getClass();
            setMethod = StsToolkit.getAccessor(objectClass, valueName, "set", valueClass);
            if(setMethod == null)
            {
                StsException.systemError("StsList.setValue() failed for valueName: " + valueName + " in class " + objectClass);
                return;
            }
            for(int n = 0; n < size; n++)
                setMethod.invoke(list[n], new Object[] { value } );
        }
        catch(Exception e)
        {
            StsException.systemError("StsList.setValue() failed for valueName: " + valueName + " in class " + objectClass);
        }
    }

    /** Returns an enumeration of the components of this StsList.
     * @return  an enumeration of the components of this list.
     */
    public final synchronized Enumeration elements()
    {
	    return new StsListEnumerator(this);
    }

    /** StsListEnumerator iterates over non-null entries in StsList */
    final class StsListEnumerator implements Enumeration
    {
        Object[] list;
        int length;
        int index;
        Object nextElement;

        StsListEnumerator(StsList objectList)
        {
            list = objectList.list;
            length = list.length;
            index = 0;
            nextElement = null;

            for(; index < length; index++)
            {
                if(list[index] != null)
                {
                    nextElement = list[index];
                    break;
                }
            }
        }

        public boolean hasMoreElements()
        {
            return nextElement != null;
        }

        public Object nextElement()
        {
            if(nextElement == null)
                return null;

            Object thisElement = nextElement;
            nextElement = null;

            for(index++; index < length; index++)
            {
                if(list[index] != null)
                {
                    nextElement = list[index];
                    break;
                }
            }
            return thisElement;
        }
    }

    /** Returns an enumeration of the components of this StsList.
     * @return  an enumeration of the components of this list.
     */
    public final synchronized Iterator iterator()
    {
	    return new StsListIterator();
    }

    /** StsListEnumerator iterates over non-null entries in StsList */
    final class StsListIterator implements Iterator
    {
        int i = 0;
        Object nextElement = null;

        StsListIterator()
        {
        }

        public boolean hasNext()
        {
            if(i >= size) return false;
            nextElement = list[i++];
            return nextElement != null;
        }

        public Object next()
        {
            return nextElement;
        }

        // this removes the last element fetched with next()
        public void remove()
        {
            if(i > 0) delete(i-1);
        }
    }

    static public Iterator emptyIterator()
    {
        StsList list = new StsList();
        return list.iterator();
    }

    public final synchronized void trimToSize()
    {
	    if (size < list.length)
        {
            Object[] newList = new Object[size];
	        System.arraycopy(list, 0, newList, 0, size);
            list = newList;
        }
	}

    public boolean listMatches(Object[] newObjects)
    {
        int nNewObjects = newObjects.length;
        if(size != nNewObjects) return false;
        for(int n = 0; n < nNewObjects; n++)
            if(list[n] != newObjects[n]) return false;
        return true;
    }

	public boolean sort()
	{
		try
		{
			if(list != null) Arrays.sort(list, 0, size);
			return true;
		}
		catch(Exception e)
		{
			StsException.systemError("Failed to sort list of class " + StsToolkit.getSimpleClassname(list[0]));
			return false;
		}
	}

    public Object[] getSortedList()
    {
        Object[] sortedList = getTrimmedList();
        Arrays.sort(sortedList, 0, size);
        return sortedList;
    }


    public StsObjectRefList convertListToRefList(StsModel model, String fieldName, StsObject parentObj)
    {
        for(int n = 0; n < size; n++)
        {
            StsObject object = (StsObject)list[n];
            if(!object.isPersistent()) object.addToModel();
        }
        return StsObjectRefList.convertListToRefList(this, fieldName, parentObj);

    }
}

