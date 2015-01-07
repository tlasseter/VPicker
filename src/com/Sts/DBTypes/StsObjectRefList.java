//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

import java.io.*;
import java.lang.reflect.*;
import java.util.*;

public class StsObjectRefList extends StsObject implements StsCustomSerializable
{
	private transient StsObjectList list;
	private transient Vector uniqueClassNames = new Vector();

	private StsObject parentObj;
	private String fieldName = "none";

	public StsObjectRefList()
	{
	}

	public StsObjectRefList(boolean persistent)
	{
		super(persistent);
	}

	private StsObjectRefList(int initSize, int incSize, String fieldName, StsObject parentObj) throws StsException
	{
		super(false);
		initialize(initSize, incSize, fieldName, parentObj, true);
	}

	private StsObjectRefList(int initSize, int incSize, String fieldName, StsObject parentObj, boolean persistent) throws StsException
	{
		super(false);
		initialize(initSize, incSize, fieldName, parentObj, persistent);
	}

    private StsObjectRefList(StsObjectList objectList, String fieldName, StsObject parentObj) throws StsException
	{
		super(false);
		this.parentObj = parentObj;
		this.fieldName = fieldName;

		if(fieldName == null)
		{
			String errorString = "StsObjectRefList.constructor() has null fieldName for " +  getParentObjectName();
			System.out.println("ERROR: " + errorString);
			throw new StsException(errorString);
		}
//		name = new String(parentObj.getName() + "[" + parentObj.index() + "]." + fieldName);

		list = objectList;
//		uniqueClassNames = new Vector();
		if (list != null)
		{
//            StsObject[] objects = list.getElements();
			int nObjects = list.getSize();
			for (int n = 0; n < nObjects; n++)
			{
				StsObject object = (StsObject) list.getElement(n);
				String className = object.getClass().getName();
				if (!uniqueClassNames.contains(className))
				{
					uniqueClassNames.add(className);
				}
			}
		}
		addToModel();
//		System.out.println("Constructed objectRefList [" + index() + "]" + " parent object: " + name);
	}

    public void initialize(int initSize, int incSize, String fieldName, StsObject parentObj, boolean persistent) throws StsException
	{
		if(fieldName == null)
		{
			String errorString = "StsObjectRefList.constructor() has null fieldName for " + getParentObjectName();
			System.out.println("ERROR: " + errorString);
			throw new StsException(errorString);
		}
		this.parentObj = parentObj;
		this.fieldName = fieldName;

		list = new StsObjectList(initSize, incSize);
		if(persistent)
        {
            addToModel();
            parentObj.fieldChanged(fieldName, this);
        }
        parentObj.objectChanged(fieldName, this);
        if(Main.isDbDebug) System.out.println("Initialized objectRefList " + toString());
	}

	private Field getField(Object obj, String fieldName) throws NoSuchFieldException
	{
		Class c = obj.getClass();
		StsDBTypeStsClass parentDBClass = (StsDBTypeStsClass)currentModel.getDatabase().getCurrentDBType(c);
		Field field = parentDBClass.getField(fieldName);
		field.setAccessible(true);
		return field;
	}

	public boolean setField(Object parentObj, String fieldName)
	{
		try
		{
			Field field = getField(parentObj, fieldName);
			if (field == null)
			{
				return false;
			}
			field.set(parentObj, this);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsObject.setField() failed for " + toString(), e, StsException.WARNING);
			return false;
		}
	}

    public Iterator getIterator()
    {
        return StsMath.arrayIterator(list.getList());
    }


    public StsObjectList.ObjectIterator getVisibleObjectIterator()
    {
        return list.getVisibleObjectIterator();
    }


    static public StsObjectRefList constructor(int initSize, int incSize, String fieldName, StsObject parentObj)
	{
		try
		{
			return new StsObjectRefList(initSize, incSize, fieldName, parentObj);
		}
		catch (Exception e)
		{
			StsException.outputException("ObjectRefList.constructor() failed for fieldname " + fieldName + " initSize=" + initSize + " incSize=" + incSize, e, StsException.WARNING);
			return null;
		}
	}

	static public StsObjectRefList constructor(int initSize, int incSize, String fieldName, StsObject parentObj, boolean persistent)
	{
		try
		{
			return new StsObjectRefList(initSize, incSize, fieldName, parentObj, persistent);
		}
		catch (Exception e)
		{
			StsException.outputException("StsObjectRefList constructor failed! initSize=" + initSize + " incSize=" + incSize, e, StsException.WARNING);
			return null;
		}
	}

    //TODO may want to consolidate StsObjectList and StsList with a flag to indicated whether list is a set of StsObjects or not 
    static public StsObjectRefList convertListToRefList(StsList list, String fieldName, StsObject parentObj)
	{
		try
		{
            StsObjectList objectList = StsObjectList.constructor(list);
            if(objectList == null) return null;
            return new StsObjectRefList(objectList, fieldName, parentObj);
		}
		catch (Exception e)
		{
			StsException.outputException("StsObjectRefList constructor failed! fieldName = " + fieldName, e, StsException.WARNING);
			return null;
		}
	}

	private StsObjectRefList(StsObject[] stsObjects, String fieldName, StsObject parentObj) throws StsException
	{
		this(new StsObjectList(stsObjects, stsObjects.length), fieldName, parentObj);
	}

	static public StsObjectRefList convertObjectsToRefList(StsObject[] stsObjects, String fieldName, StsObject parentObj)
	{
		try
		{
			return new StsObjectRefList(stsObjects, fieldName, parentObj);
		}
		catch (Exception e)
		{
			StsException.outputException("StsObjectRefList constructor failed! fieldName = " + fieldName, e, StsException.WARNING);
			return null;
		}
	}
	public boolean initialize(StsModel model)
	{
		return true;
	}

	public boolean add(StsObject obj)
	{
		if (obj == null)
		{
			return false;
		}
		if (getIndex() >= 0 && obj.getIndex() < 0)
		{
			StsException.systemError("Can't add non-persistent object " + obj.getLabel() + " to StsObjectRefList owned by  " + parentObj.toString());
			return false;
		}
		int listIndex = list.add(obj);
		if (listIndex < 0) return false; // failed to add object

		if(getIndex() < 0) return true; // don't create db command below for nonpersistent list

		StsAddRefCmd cmd = new StsAddRefCmd(this, obj, listIndex);
		currentModel.addTransactionCmd("addRefListObj", cmd);
//        StsTransaction.add(cmd);

		currentModel.getDatabase().getCurrentDBType(obj.getClass());
		String className = obj.getClass().getName();
		if (!uniqueClassNames.contains(className))
		{
			uniqueClassNames.addElement(className);
		}
		return true;
	}

	public boolean addDB(StsObject obj)
	{
		if (obj == null) return false;
		int listIndex = list.add(obj);
		return listIndex >= 0;
	}

	public void add(StsObject[] objs)
	{
		if (objs == null)
		{
			return;
		}
		for (int i = 0; i < objs.length; i++)
		{
			add(objs[i]);
		}
	}

	public boolean deleteElement(Object obj)
	{
		return list.delete(obj) >= 0;
	}

	private void unInsert(StsObject obj) throws StsException
	{
		list.unInsert(obj);
	}

	public boolean delete(StsObject obj)
	{
		if (obj == null || obj.getIndex() < 0)return false;
        int listIndex = list.getIndex(obj);
		if(list.delete(obj) == StsParameters.NO_MATCH) return false;
        StsDeleteRefCmd cmd = new StsDeleteRefCmd(this, obj, listIndex);
        currentModel.addTransactionCmd("deleteRefListObj", cmd);
		return true;
	}

	/** delete object at this sequential location */
	public boolean delete(int index)
	{
		return list.delete(index);
	}

	/** delete object which has this index */
	public boolean deleteObjectWithIndex(int index)
	{
		int size = list.getSize();
		for (int n = 0; n < size; n++)
		{
			StsObject object = (StsObject) list.getElement(n);
			if (object.getIndex() == index)
			{
				list.delete(n);
				return true;
			}
		}
		return false;
	}

	public StsObject getElementWithIndex(int index)
	{
		int size = list.getSize();
		for (int n = 0; n < size; n++)
		{
			StsObject object = (StsObject) list.getElement(n);
			if (object.getIndex() == index)return object;
		}
		return null;
	}

	/** Delete all elements in this list and deletes this list.
	 * Reverse order avoids some reindexing.
	 */
	public void deleteAll()
	{
        // delete each list entry
		int nObjs = getSize();
		for (int n = nObjs - 1; n >= 0; n--)
		{
			StsObject obj = (StsObject) list.getElement(n);
			if(obj != null) obj.delete();
		}
        list.deleteAll(); // removes all entries from list
		uniqueClassNames.clear();
		delete(); /** delete this list */
	}

	/** delete all objects in the list from the db,
	 *  and clear objects from list, but keep list.
	 *  Reverse order avoids some reindexing.
	 */
	public void deleteAllElements()
	{
		int nObjects = getSize();
		for (int n = 0; n < nObjects; n++)
		{
			StsObject obj = (StsObject) list.getElement(n);
			obj.delete();
		}
		list.deleteAll();
	}

	static public void deleteAll(StsObjectRefList refList)
	{
		if (refList != null)
		{
			refList.deleteAll();
		}
	}

	public void insertBefore(int index, StsObject newObj)
	{
		int insertedIndex = list.insertBefore(index, newObj);
		// Entire list has been changed: need Cmd and Transaction here
		StsAddRefCmd cmd = new StsAddRefCmd(this, newObj, insertedIndex);
		currentModel.addTransactionCmd("insertBeforeObj", cmd);
	}

	public void insertAfter(int index, StsObject newObj)
	{
		int insertedIndex = list.insertAfter(index, newObj);
		// Entire list has been changed: need Cmd and Transaction here
		StsAddRefCmd cmd = new StsAddRefCmd(this, newObj, insertedIndex);
		currentModel.addTransactionCmd("insertAfterObj", cmd);
	}

	public void insertBefore(StsObject obj, StsObject newObj)
	{
		int index = list.getIndex(obj);
		insertBefore(index, newObj);
	}

	public void insertAfter(StsObject obj, StsObject newObj)
	{
		int index = list.getIndex(obj);
		insertAfter(index, newObj);
	}

	public int getIndex(StsObject obj)
	{
		return list.getIndex(obj);
	}

	public boolean insertObject(StsObject obj, int index)
	{
		int size = getSize();
		if (index > size)
		{
            list.resize(index+1);
            list.setElement(obj, index);
//            String className = obj.getClass().getName();
//			StsException.systemError("          insert object failed for " + className + "[" + obj.getIndex() + "]. Current list size is " + size);
//			return false;
		}
		else if (index == size)
		{
			add(obj);
		}
		else if (index < size)
		{
			StsObject currentObj = (StsObject) list.getElement(index);
			if (currentObj == obj || currentObj == null)
			{
				list.setElement(obj, index);
				if (Main.isDbDebug)
				{
					String className = obj.getClass().getName();
					if (currentObj == null)
					{
						System.out.println("          replacing null object at index " + index + " with " + className + "[" + obj.getIndex() + "]");
					}
					else
					{
						System.out.println("          replacing existing object at index " + index + " with " + className + "[" + obj.getIndex() + "]");
					}
				}
			}
			else
			{
				list.insertBefore(index, obj);
				if (Main.isDbDebug)
				{
					String className = obj.getClass().getName();
					System.out.println("          inserting before object at index " + index + " with" + className + "[" + obj.getIndex() + "]");
				}
			}
			return true;
		}
		return true;
	}

	public StsObject setCurrentPosition(int index)
	{
		return (StsObject) list.setCurrentPosition(index);
	}

	public boolean hasObject(StsObject obj)
	{
		return list.contains(obj);
	}

	public boolean hasObjectAtIndex(StsObject obj, int index)
	{
		if (index >= getSize()) return false;
		return obj == list.getElement(index);
	}

    public StsObject getObjectWithName(StsObject otherObject)
    {
        if(otherObject == null) return null;
        return list.getObjectWithName(otherObject.getName());
    }

	public int getSize()
	{
		if (list == null)return 0;
		return list.getSize();
	}

	public StsObject getElement(int index)
	{
		return (StsObject) list.getElement(index);
	}

	public StsObject getFirst()
	{
		return (StsObject) list.getFirst();
	}

	public StsObject getNext()
	{
		return (StsObject) list.getNext();
	}

	public StsObject getLast()
	{
		return (StsObject) list.getLast();
	}

	public StsObject getSecondToLast()
	{
		return (StsObject) list.getSecondToLast();
	}

	public StsObject getSecond()
	{
		return (StsObject) list.getSecond();
	}

	public Enumeration elements()
	{
		return list.elements();
	}

	public boolean contains(StsObject object)
	{
		return list.contains(object);
	}

	public StsObjectList getList()
	{
		return list;
	}

    public void clear()
    {
        getList().clear();
    }

    public Object[] getTrimmedList()
    {
        return list.getTrimmedList();
    }

	public Object[] getArrayList()
	{
		return list.getList();
	}

	public StsObject[] getListCopy()
	{
		return list.getListCopy();
	}

	public Object getCastListCopy()
	{
		return list.getCastListCopy();
	}

	static public int getSize(StsObjectRefList list)
	{
		if (list == null)
		{
			return 0;
		}
		else
		{
			return list.getSize();
		}
	}

	public Object copyArrayList()
	{
		return list.copyArrayList();
	}

	public Object copyArrayList(Class componentType)
	{
		return list.copyArrayList(componentType);
	}

    public StsObjectRefList deepCopy()
    {
        if(list == null) return null;
        StsObjectList newList = list.deepCopy();
        StsObjectRefList newObjectRefList = (StsObjectRefList)clone();
        newObjectRefList.list = newList;
        return newObjectRefList;
    }

    // handles null objects
	static public int staticGetSize(StsObjectRefList objectRefList)
	{
		if (objectRefList == null)
		{
			return 0;
		}
		return objectRefList.getSize();
	}

	public void forEach(String methodName)
	{
		list.forEach(methodName, null);
	}

	/** If a single argument, make an array of length 1 for compatibility */
	public void forEach(String methodName, Object arg)
	{
		list.forEach(methodName, new Object[] {arg});
	}

	public void forEach(String methodName, Object[] args)
	{
		list.forEach(methodName, args);
	}

	/** return an array of names */
	public String[] getNames()
	{
		if (list == null)
		{
			return null;
		}
		return list.getNames();
	}

	/** return an object from its name */
	public StsObject getElement(String name)
	{
		if (list == null)
		{
			return null;
		}
		return list.getElement(name);
	}

	/** return an array of objects */
	public StsObject[] getElements()
	{
		return list.getElements();
	}

    public boolean listMatches(StsObject[] newObjects)
    {
        return list.listMatches(newObjects);
    }

	/** returns an array cast to class */
	public Object getCastList()
	{
		return StsMath.arraycopy(list.getList(), list.getSize());
	}

	public Object getCastList(Class componentClass)
	{
		Object[] objects = list.getList();
		return StsMath.arraycopy(objects, getSize(), componentClass);
	}

	public void reverseOrder()
	{
		list.reverseOrder();
	}

	public boolean sort()
	{
		return list.sort();
	}

    public Object[] getSortedList()
    {
        return list.getSortedList();
    }

    public int getParentIndex() { return parentObj.getIndex(); }
	public StsObject getParentObject() { return parentObj; }
	public String getParentObjectName() { return parentObj.getClass().getName() + "[" + parentObj.getIndex() + "]"; }

	/** return name of this list */
	public String getName()
	{
        return toString();
	}

    public String toString()
    {
        if(this.parentObj != null)
            return StsToolkit.getSimpleClassname(this) + " " + getParentObjectName() + "." + fieldName;
        else
            return StsToolkit.getSimpleClassname(this) + " " + "NullParent" + "." + fieldName;
    }

	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException
	{
		int numClasses = uniqueClassNames.size();
		int size = list.getSize();
		out.writeInt(size);
		if (size > 0)
		{
			out.writeInt(list.getIncSize());
			out.writeBoolean(numClasses == 1);
			// write out the index of all the non-null entries
			if (numClasses == 1) // homogeneous list
			{
				StsDBTypeClass dbType = (StsDBTypeClass)out.getOutputDBType(list.getElement(0).getClass());
				out.writeInt(dbType.getIndex());
				for (int n = 0; n < size; n++)
				{
					StsObject obj = (StsObject)list.getElement(n);
					out.writeInt(obj.getIndex());
				}
			}
			// write out the classIndex and index of all the non-null entries
			else
			{
				for (int n = 0; n < size; n++)
				{
					StsObject obj = (StsObject)list.getElement(n);
					StsDBTypeClass dbType = (StsDBTypeClass)out.getOutputDBType(obj.getClass());
					out.writeInt(dbType.getIndex());
					out.writeInt(obj.getIndex());
				}
			}
		}
	}

    public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException
	{
		int size = in.readInt();
		if (size < 1)
		{
			list = new StsObjectList(1);
		}
		else
		{
			int listIncSize = in.readInt();
			// read/resolve the references
			list = new StsObjectList(size, listIncSize);
			list.setSize(size);
			if (in.readBoolean()) // homogeneous list
			{
				int dbClassTypeIndex = in.readInt();
				StsDBTypeClass dbType = (StsDBTypeClass)in.getInputDBType(dbClassTypeIndex);
				if (dbType == null)
				{
					System.out.println("Instancelist not found for " + dbClassTypeIndex);
				}
				else
				{
					for (int i = 0; i < size; i++)
					{
						Object obj = in.getModelObject(dbType.getClassType(), in.readInt());
						list.list[i] = obj;
						String className = obj.getClass().getName();
						if (!uniqueClassNames.contains(className))
						{
							uniqueClassNames.addElement(className);
						}
					}
				}
			}
			else // heterogeneous list
			{
				for (int i = 0; i < size; i++)
				{
					int classIndex = in.readInt();
					int objIndex = in.readInt();
					StsDBTypeClass dbType = (StsDBTypeClass)in.getInputDBType(classIndex);
					Object obj = in.getModelObject(dbType.getClassType(), objIndex);
                    if(obj == null)
                    {
                        StsException.systemError(this, "readObject", "obj is null.");
                        continue;
                    }
                    list.list[i] = obj;
					String className = obj.getClass().getName();
					if (!uniqueClassNames.contains(className))
					{
						uniqueClassNames.addElement(className);
					}
				}
			}
		}
	}

    public void exportObject(StsDBFileObjectTrader objectTrader) throws IllegalAccessException
	{
        int size = list.getSize();
        for (int n = 0; n < size; n++)
        {
            StsObject obj = (StsObject)list.getElement(n);
            StsDBTypeObject dbType = (StsDBTypeObject)objectTrader.getCurrentDBType(obj);
            dbType.exportObject(objectTrader, obj);
        }
	}
}
