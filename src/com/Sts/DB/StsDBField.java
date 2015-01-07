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

public class StsDBField implements Comparable, StsSerializable
{
	/** name of this field in class */
	private String fieldName;
	/** type of this field as string; e.g., int[][] or java.lang.Object */
	private String typeName;
	/** Object type of this field */
	private transient StsDBType dbType = null;
	/** set from the class description */
	private transient Field field;
	/** if true, clone oldObject for this field: simply set reference to oldObject  */
	private transient boolean cloneStsObject = false;

	static final protected boolean debug = Main.isDbIODebug;

	public StsDBField()
	{
		super();
	}

	private StsDBField(Field field)
	{
		this.field = field;
		field.setAccessible(true);
		Class fieldClass = field.getType();
		typeName = fieldClass.getName();
		fieldName = field.getName();
	}

	static public StsDBField constructor(Field field)
	{
		try
		{
			int mod = field.getModifiers();
			if (Modifier.isStatic(mod) || Modifier.isTransient(mod) || field.getName().equals("this$0"))
				return null;
			return new StsDBField(field);
		}
		catch(Exception e)
		{
			StsException.outputException("DBFieldType.constructor() failed for field " + field.getName(), e, StsException.WARNING);
			return null;
		}
	}

	public Field getField()
	{
		return field;
	}

	public String getFieldName()
	{
		return fieldName;
	}

	public StsDBType getDBType()
	{
		return dbType;
	}

	/** Allow the field to be set from another field. This is to allow template classes to read
	 *  objects from disk.
	 */
	protected void setField(StsDBField dbField)
	{
		this.field = dbField.field;
	}

	public String getFieldClassName()
	{
		return typeName;
	}

	/** Get the value of this field in this object and write it out */
	public void write(StsDBOutputStream out, Object classObject) throws IOException, IllegalAccessException
	{
        try
        {
            Class fieldType = field.getType();
            Object fieldObject = field.get(classObject);
            if (fieldObject == null || fieldType.isPrimitive() || fieldType == String.class)
            {
                dbType.writeField(out, classObject, field);
            }
            else
            {
                StsDBType actualDBType = out.getOutputDBType(fieldObject);
                actualDBType.writeField(out, classObject, field);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "write(out, classObject)", "class: " + StsToolkit.getSimpleClassname(classObject) + " fieldName: " + fieldName, e);
        }
    }

	/** Get the value of this field in this object and write it out */
	public void export(StsDBFileObjectTrader objectTrader, Object classObject) throws IllegalAccessException
	{
		Object o = field.get(classObject);
		if (o == null || field.getType().isPrimitive() || field.getType() == String.class)
		{
			return;
		}
		else
		{
			StsDBType actualDBType = objectTrader.getCurrentDBType(field.get(classObject));
			actualDBType.exportField(objectTrader, classObject, field);
		}
	}

	public void read(StsDBInputStream in, Object classObject) throws IllegalAccessException, IOException
	{
        if(field == null)
        {
            StsException.systemError(this, "read", "No current field exists in class " +
                StsToolkit.getSimpleClassname(classObject) + " for field name: " + fieldName);
        }
		dbType.readField(in, classObject, field);
	}

	public int compareTo(Object otherField)
	{
		String otherTypeName = ((StsDBField)otherField).typeName;
		if(typeName.equals(otherTypeName))
			return 0;
		else
			return -1;
	}

	public void copy(Object oldObject, Object newObject)  throws IOException, IllegalAccessException
	{
		Object oldFieldObject = field.get(oldObject);
		if(oldFieldObject == null) return;
		if(dbType instanceof StsDBTypeClass)
		{
			StsDBType fieldObjectType = getDBType(oldFieldObject, getCopyTypeList());
			fieldObjectType.copyField(oldObject, newObject, field);
		}
		else
			dbType.copyField(oldObject, newObject, field);
	}

	static private StsDBObjectTypeList getCopyTypeList()
	{
		return StsObjectCopier.typeList;
	}

	private StsDBType getDBType(Object object, StsDBObjectTypeList typeList)
	{
		return typeList.getDBType(object);
	}

	public void clone(Object oldObject, Object newObject) throws IllegalAccessException
	{
		if(debug) System.out.println("cloning field " + fieldName + " for class " + typeName);
		((StsDBTypeClass)dbType).cloneField(oldObject, newObject, field);
	}

	public void addToModel(Object object, StsModel model)
	{
		Object fieldObject = null;

		if (cloneStsObject)
			return;
		try
		{
			fieldObject = field.get(object);
			if (fieldObject == null) return;
			dbType.addToModel(fieldObject, model);
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "addToModel", "Failed for object " + StsToolkit.getSimpleClassname(fieldObject), e);
		}
	}

	public boolean equals(StsDBField other)
	{
		return this.typeName.equals(other.typeName);
	}

	public boolean initializeFieldTypes(StsDBObjectTypeList typeList, Class parentClassType)
	{
		dbType = typeList.getDBType(typeName, parentClassType);
		return dbType != null;
	}

	public void setCloneFlag(boolean clone)
	{
		if(debug) System.out.println("fieldname clone flag set to " + clone);
		cloneStsObject = clone;
	}

	public boolean getCloneStsObject()
	{
		return cloneStsObject;
	}

    public String toString() { return fieldName + " " + typeName; }
}
