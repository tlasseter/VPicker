package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsDBObjectTypeList
{
	private Map dbTypesByName = new TreeMap();
	private StsDBType[] dbTypesByIndex = new StsDBType[200];
	private int classIndex = 0;
	private int arrayTypeIndex = -1;
	private static int listIncSize = 50;
    private String name = "none";
//    private boolean dbTypeTemporary = false;

    private ArrayList newDBTypeClasses = new ArrayList();

    private boolean debug = Main.isDbCmdDebug;

    public StsDBObjectTypeList(String name)
	{
        this.setName(name);
        arrayTypeIndex = getNextClassIndex();
		store(new StsDBTypeArray(arrayTypeIndex));

		store(new StsDBTypeInt(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Integer.class, new StsDBObjectIOInteger()));

		store(new StsDBTypeFloat(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Float.class, new StsDBObjectIOFloat()));

		store(new StsDBTypeDouble(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Double.class, new StsDBObjectIODouble()));

		store(new StsDBTypeBoolean(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Boolean.class, new StsDBObjectIOBoolean()));

		store(new StsDBTypeByte(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Byte.class, new StsDBObjectIOByte()));

		store(new StsDBTypeLong(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Long.class, new StsDBObjectIOLong()));

		store(new StsDBTypeShort(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Short.class, new StsDBObjectIOShort()));

		store(new StsDBTypeChar(getNextClassIndex()));
		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Character.class, new StsDBObjectIOCharacter()));

		store(new StsDBTypePrimitiveObject(getNextClassIndex(), String.class, new StsDBObjectIOString()));

		store(new StsDBTypePrimitiveObject(getNextClassIndex(), Class.class, new StsDBObjectIOClass()));

		store(new StsDBTypeClass(getNextClassIndex(), Object.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBField.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBTypeClass.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBTypeStsSerializable.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBTypeArray.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBTypeStsClass.class));
		store(new StsDBTypeStsSerializable(getNextClassIndex(), StsDBTypeModel.class));

		storeJavaClass(new StsDBTypeJavaColorSpace(getNextClassIndex()));
		storeJavaClass(new StsDBTypeJavaColor(getNextClassIndex()));
		storeJavaClass(new StsDBTypeJavaLocale(getNextClassIndex()));

		storeJavaClass(new StsDBTypeClass(getNextClassIndex(), java.awt.Point.class));
		storeJavaClass(new StsDBTypeClass(getNextClassIndex(), java.awt.Dimension.class));
		storeJavaClass(new StsDBTypeClass(getNextClassIndex(), java.awt.Rectangle.class));

//		classes = new LinkedList();
	}

	public StsDBObjectTypeList(StsDBObjectTypeList inputList)
	{
		this("inputCopy");
        Iterator itr = inputList.dbTypesByName.values().iterator();
//        Iterator itr = inputList.dbClassesByName..iterator();
		while (itr.hasNext())
		{
			StsDBType dbType = (StsDBType)itr.next();
            if (dbType != null && !dbType.isTemporary())
            {
                StsDBType currentDBType = getDBType(dbType.getClassType());
                if(currentDBType != null)
                {
                    if(currentDBType.getIndex() == dbType.getIndex()) continue;
                    move(currentDBType, dbType.getIndex());
                    currentDBType.setMatchesOnDiskVersion(true);
                }
                else
                {
                    insert(dbType, dbType.getIndex());
				    dbType.setMatchesOnDiskVersion(true);
                }
            }
		}
		classIndex = getFirstFreeIndex();
	}

    public void removeTemporaryTypes()
    {
        int nClasses = dbTypesByIndex.length;
        StsDBType[] newDbTypesByIndex = new StsDBType[nClasses];
        for(int n = 0; n < dbTypesByIndex.length; n++)
        {
            if(dbTypesByIndex[n] != null)
            {
                StsDBType dbType = dbTypesByIndex[n];
                if(dbType == null) continue;
                if(dbType.isTemporary())
                {
                    dbTypesByIndex[n] = null;
                    dbTypesByName.remove(dbType.getTypeName());
                }
                else
                    newDbTypesByIndex[n] = dbType;
            }
        }
        dbTypesByIndex = newDbTypesByIndex;
        classIndex = getFirstFreeIndex();
        clearNewClasses();
    }

	public StsDBType getDBType(Class c)
	{
        if(c == null) return null;
        try
        {
            String className = c.getName();
            StsDBType dbType = (StsDBType) dbTypesByName.get(className);
            if (dbType == null)
            {
                dbType = addDBType(c, null);
            }
            return dbType;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getDBType", e);
            return null;
        }
    }

	public StsDBType getDBType(Object obj)
	{
		Class c = obj.getClass();
		return getDBType(c);
	}

	public StsDBType getDBTypeForceInsert(StsDBTypeClass dbTypeClass)
	{
		StsDBType dbType = (StsDBType) dbTypesByName.get(dbTypeClass.getTypeName());
		if (dbType == null)
		{
			int index = dbTypeClass.getIndex();
			dbType = addDBType(dbTypeClass.getClassType(), index);
		}
		else
		{
			move(dbType, dbTypeClass.getIndex());
		}
		return dbType;
	}

	public StsDBType getDBType(String className, Class parentClassType)
	{
		StsDBType dbClass = (StsDBType) dbTypesByName.get(className);
		if (dbClass == null)
		{
			try
			{
				dbClass = addDBType(Class.forName(className), parentClassType);
                if(!dbClass.isArray() && !newDBTypeClasses.contains(dbClass))
        		    newDBTypeClasses.add(dbClass);
            }
			catch (ClassNotFoundException ex)
			{
				StsException.systemError(this, "getDBType(String, Class) failed trying to find class " + className);
			}
		}
		return dbClass;
	}

	public StsDBType getDBType(int index)
	{
		if (index == this.arrayTypeIndex)
		{
			return StsDBTypeArray.getArrayTypeInstance(arrayTypeIndex);
		}
		StsDBType dbt = dbTypesByIndex[index];
		return dbt;
	}
/*
	public Object[] getNonDefaultDBClassArray()
	{
		return classes.toArray();
	}
*/
	private void store(StsDBType dbType)
	{
		int index = dbType.getIndex();
		dbTypesByName.put(dbType.getTypeName(), dbType);
		StsMath.arrayInsertElement(dbTypesByIndex, dbType, index, listIncSize);
		classIndex = getFirstFreeIndex();
		dbType.initializeFieldTypes(this);
	}

    private void store(StsDBTypeClass dbClassType)
	{
		int index = dbClassType.getIndex();
		dbTypesByName.put(dbClassType.getTypeName(), dbClassType);
		StsMath.arrayInsertElement(dbTypesByIndex, dbClassType, index, listIncSize);
		classIndex = getFirstFreeIndex();
		dbClassType.initializeFieldTypes(this);
		if(!dbClassType.isArray()) newDBTypeClasses.add(dbClassType);
        if(debug) System.out.println("Stored dbClassType " + dbClassType.getTypeName() + " in list " + name);
    }

    private void storeJavaClass(StsDBTypeClass dbClassType)
	{
		int index = dbClassType.getIndex();
		dbTypesByName.put(dbClassType.getTypeName(), dbClassType);
		StsMath.arrayInsertElement(dbTypesByIndex, dbClassType, index, listIncSize);
		classIndex = getFirstFreeIndex();
		dbClassType.initializeFieldTypes(this);
        if(debug) System.out.println("Stored JavaClassType " + dbClassType.getTypeName() + " to list " + name);
    }

    private void move(StsDBType oldDBClassType, int newIndex)
	{
		int oldIndex = oldDBClassType.getIndex();
		if (oldIndex == newIndex)
		{
			return;
		}
		dbTypesByIndex[oldIndex] = null;
		if (dbTypesByIndex.length > newIndex && dbTypesByIndex[newIndex] != null)
		{
			dbTypesByName.remove(dbTypesByIndex[newIndex].getTypeName());
			dbTypesByIndex[newIndex] = null;
		}
        oldDBClassType.setIndex(newIndex);
        StsMath.arrayInsertElement(dbTypesByIndex, oldDBClassType, newIndex, listIncSize);
		newDBTypeClasses.add(oldDBClassType);
        if(debug) System.out.println("Moved dbClassType " + oldDBClassType.getTypeName() + " from index " + oldIndex + " to " + newIndex + " in list " + name);
    }

	protected void insert(StsDBType dbClassType, int newIndex)
	{
		int index = dbClassType.getIndex();
		StsDBType existingDBType = (StsDBType) dbTypesByName.get(dbClassType.getTypeName());
		if (existingDBType != null)
		{
			dbTypesByIndex[existingDBType.getIndex()] = null;
			dbTypesByName.remove(existingDBType.getTypeName());
			newDBTypeClasses.remove(existingDBType);
		}
		if (dbTypesByIndex.length > index && dbTypesByIndex[index] != null)
		{
			dbTypesByName.remove(dbTypesByIndex[index].getTypeName());
			dbTypesByIndex[index] = null;
		}
		dbTypesByName.put(dbClassType.getTypeName(), dbClassType);
		StsMath.arrayInsertElement(dbTypesByIndex, dbClassType, newIndex, listIncSize);
		dbClassType.setIndex(newIndex);
		classIndex = getFirstFreeIndex();
		dbClassType.initializeFieldTypes(this);
		newDBTypeClasses.add(dbClassType);
         if(debug) System.out.println("Inserted dbClassType " + dbClassType.getTypeName() + " in list " + name);
    }

    public void clearNewClasses()
    {
        newDBTypeClasses.clear();
    }

    public Object[] getNewClasses() { return newDBTypeClasses.toArray(); }

    public void debugPrintNewClasses()
    {
        Iterator iter = newDBTypeClasses.iterator();
        System.out.println("new dbTypeClasses");
        while(iter.hasNext())
        {
            StsDBType dbType = (StsDBType)iter.next();
            System.out.println("    " + dbType.getTypeName());
        }
    }

    public void printDBTypes()
    {
        System.out.println("dbTypeByIndex");
        for(int n = 0; n < dbTypesByIndex.length; n++)
        {
            StsDBType dbType = dbTypesByIndex[n];
            if(dbType == null) continue;
            System.out.println("    " + dbType.getTypeName());
        }
    }

    private int getFirstFreeIndex()
	{
		int result = dbTypesByIndex.length;
		while (dbTypesByIndex[result - 1] == null)
		{
			result--;
		}
		return result;
	}

	private StsDBType addDBType(Class c, Class parentClassType)
	{
		if (c.isArray())
			return StsDBTypeArray.getArrayTypeInstance(arrayTypeIndex, c, this);
		StsDBTypeClass dbClass = constructor(c, getNextClassIndex(), parentClassType);
		store(dbClass);
		return dbClass;
	}

	private StsDBType addDBType(Class c, int newIndex)
	{
        if(c == null) return null;
        if (c.isArray())
			return StsDBTypeArray.getArrayTypeInstance(arrayTypeIndex, c, this);
		StsDBType dbClass = constructor(c, newIndex, null);
		insert(dbClass, newIndex);
		return dbClass;
	}


    private StsDBType checkAddDBType(Class c, int newIndex)
    {
        if(c == null) return null;

        if (c.isArray())
            return StsDBTypeArray.getArrayTypeInstance(arrayTypeIndex, c, this);
        StsDBType dbClass = constructor(c, newIndex, null);
        insert(dbClass, newIndex);
        return dbClass;
    }


    private StsDBTypeClass constructor(Class fieldClass, int newIndex, Class parentClassType)
	{
		if (fieldClass.isArray())
		{
			System.out.println("Should never get here.");
		}
		if (StsModel.class.isAssignableFrom(fieldClass))
			return new StsDBTypeModel(newIndex);
		if (StsObject.class.isAssignableFrom(fieldClass))
			return new StsDBTypeStsClass(newIndex, fieldClass);
		if (StsSerializable.class.isAssignableFrom(fieldClass))
			return new StsDBTypeStsSerializable(newIndex, fieldClass);
		if (fieldClass.getPackage().getName().startsWith("java"))
		{
			StsException.systemError(this, "constructor(fieldClass, newIndex, parentClass)",
                    "Warning! Adding java package dbType for field " + fieldClass.getName() + " for class: " + parentClassType +
                    ". Probably should make this field transient.");
		}
		return new StsDBTypeClass(newIndex, fieldClass);
	}

	private synchronized int getNextClassIndex()
	{
		return classIndex++;
	}

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }
}
