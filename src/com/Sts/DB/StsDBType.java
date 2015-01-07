package com.Sts.DB;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;

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

public abstract class StsDBType implements StsSerializable
{
	static final public boolean debug = Main.isDbIODebug;
	private int index = -99;
	private String typeName;
    protected boolean isTemporary = false;
    protected transient Class classType;   

    public StsDBType()
	{
	}

	public StsDBType(Class classType)
	{
		super();
		this.classType = classType;
		typeName = classType.getName();
    }

	public StsDBType(int index, Class classType)
	{
		super();
		this.index = index;
		this.classType = classType;
		typeName = classType.getName();
    }

	public StsDBType(int index)
	{
		super();
		this.index = index;
		this.classType = getClass();
		typeName = classType.getName();
    }
/*
    static public void setTypeIsTemporary(boolean temporary) { typeIsTemporary = temporary; }

    public boolean isTemporary() { return isTemporary; }
*/
    protected void initializeAfterLoad(StsDBObjectTypeList typeList) throws ClassNotFoundException
	{
		try
		{
            if(typeName == null) return;
            classType = Class.forName(typeName);
		}
		catch (Exception ex)
		{
			StsException.systemError(this, "initializeAfterLoad", "Failed to find class " + typeName);
		}
	}

	protected void initializeFieldTypes(StsDBObjectTypeList typeList)
	{
	}

	protected void setClassFromTypeName(String typeName) throws ClassNotFoundException
	{
		classType = Class.forName(typeName);
		this.typeName = typeName;
	}

	protected void setClassFromObject(Object object)
	{
		classType = object.getClass();
		this.typeName = classType.getName();
	}

	abstract public void writeField(StsDBOutputStream out, Object classObject, Field field) throws IllegalAccessException, IOException;
	abstract public void readField(StsDBInputStream in, Object classObject, Field field) throws IllegalAccessException, IOException;
	abstract public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException, IllegalAccessException;
	abstract public Object readArrayValues(StsDBInputStream in, Object arrayObject) throws IOException;

	abstract public void copyField(Object oldObject, Object newObject, Field field) throws IOException, IllegalAccessException;
	abstract public Object copyArrayValues(Object arrayObject) throws IOException, IllegalAccessException;

	public void exportField(StsDBFileObjectTrader objectTrader, Object obj, Field field) throws IllegalAccessException
	{
	}

	public void exportArrayValues(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
	}

	/** DBClassType returns false. */
	public boolean isDBClass()
	{
		return false;
	}

    public boolean isArray() { return false; }

    public int getIndex()
	{
		return index;
	}

	public void setIndex(int index)
	{
		this.index = index;
	}

	public Class getClassType()
	{
		if (classType != null) return classType;
        if(typeName == null) return null;
        try
        {
            classType = Class.forName(typeName);
        }
        catch (ClassNotFoundException ex)
        {
            classType = null;
        }
		return classType;
	}

	public String getTypeName()
	{
		return typeName;
	}

	public boolean getMatchesOnDiskVersion()
	{
		return true;
	}

	public void setMatchesOnDiskVersion(boolean matches)
	{
	}

    public void setMatchesOnDiskVersion(StsDBType onDiskVeriosn)
	{
	}

	/** overridden in subclasses which handle StsObjects (DBStsClassType) and StsObject arrays (DBArrayType) */
	public void addToModel(Object object, StsModel model)
	{
	}

	public boolean representsStsSerialisableObject()
	{
		return false;
	}

	public boolean equals(StsDBType other)
	{
		return this.typeName.equals(other.typeName);
	}
    public boolean isTemporary() { return isTemporary; }
}
