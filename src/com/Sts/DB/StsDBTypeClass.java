package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;
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

public class StsDBTypeClass extends StsDBTypeObject implements StsSerializable
{
	static public final int NULL_REFERENCE = -99;
	/** fields for this class.  TreeMap is not directly serialized; converted to Object array and then serialized. */
	private StsDBField[] dbFields;
	/** fields for this class.  TreeMap is not directly serialized; converted to array and then serialized. */
	private transient TreeMap dbFieldsMap;
	private transient boolean matchesOnDiskVersion = false;
    {
        isTemporary = typeIsTemporary;
    }

    static boolean typeIsTemporary = false;
    public StsDBTypeClass()
	{
		super();
	}

	public StsDBTypeClass(int index, Class classType)
	{
		super(index, classType);
		dbFieldsMap = getAllFields();
		setFieldsFromTreeMap();
        if(isTemporary)
        {
            if(Main.isDbIODebug)
            {
                System.out.println("----Created temporary dbType: " + getTypeName());
                if(getTypeName().equals(StsWell.class.getName()))
                    System.out.println("Temp well");
            }
        }
	}

    public StsDBTypeClass(int index)
	{
		super(index);
		dbFieldsMap = getAllFields();
		setFieldsFromTreeMap();
        if(isTemporary)
        {
            if(Main.isDbIODebug)
            {
                System.out.println("----Created temporary dbType: " + getTypeName());
                if(getTypeName().equals(StsWell.class.getName()))
                    System.out.println("Temp well");
            }
        }
	}

    public StsDBField getDBField(String fieldName)
	{
		StsDBField dbField = (StsDBField)dbFieldsMap.get(fieldName);
		if (dbField == null)
		{
			return null;
		}
		return dbField;
	}

	public Field getField(String fieldName)
	{
		StsDBField dbField = (StsDBField)dbFieldsMap.get(fieldName);
		if (dbField == null)
		{
			return null;
		}
		return dbField.getField();
	}

	public Field getField(int fieldNumber)
	{
        if(dbFields == null) return null;
        if(fieldNumber < 0 || fieldNumber > dbFields.length-1)
        {
            StsException.systemError(this, "getField", " Number of fields: " + dbFields.length + ".  Requested fieldNumber " + fieldNumber);
            return null;
        }
        return dbFields[fieldNumber].getField();
	}

	public int getIndexOfField(String fieldName)
	{
		StsDBField dbField = (StsDBField)dbFieldsMap.get(fieldName);
		if (dbField == null)
		{
			return -1;
		}

		for (int i = 0; i < dbFields.length; i++)
		{
			StsDBField thisField = dbFields[i];
			if (dbField == thisField)
			{
				return i;
			}
		}
		return -1;
	}

	public StsDBField getDBFieldType(String fieldName)
	{
		return (StsDBField)dbFieldsMap.get(fieldName);
	}

	public void setFieldsFromCurrentDBClass(StsDBTypeClass dbClass)
	{
		// 'this' is the template dbClass (on disk definition).
		// 'that' is the class representing the actual object.

		setFieldsFromTreeMap();
        if(dbFields == null) return;
        for (int i = 0; i < dbFields.length; i++)
		{
			StsDBField thisField = (StsDBField)dbFields[i];
			StsDBField thatField = (StsDBField)dbClass.dbFieldsMap.get(thisField.getFieldName());
			if(thatField != null && thisField.equals(thatField))
			{
				thisField.setField(thatField);
			}
		}

		dbClass.matchesOnDiskVersion = same(dbClass);
	}

	private boolean same(StsDBTypeClass that)
	{
        if(that.dbFields == null)
        {
            return this.dbFields == null || this.dbFields.length == 0;
        }
        if (this.dbFields.length != that.dbFields.length)
		{
			return false;
		}
		for (int i = 0; i < this.dbFields.length; i++)
		{
			StsDBField thisField = dbFields[i];
			StsDBField thatField = (StsDBField)that.dbFieldsMap.get(thisField.getFieldName());
			if (thatField == null)
			{
				return false;
			}
			if (! thisField.getFieldClassName().equals(thatField.getFieldClassName()))
			{
				return false;
			}
		}
	   return true;
	}

	public void setMatchesOnDiskVersion(StsDBType onDiskVeriosn)
	{
		matchesOnDiskVersion = same((StsDBTypeClass)onDiskVeriosn);
	}

	public void setMatchesOnDiskVersion(boolean matches)
	{
		matchesOnDiskVersion = matches;
	}

    static public void setTypeIsTemporary(boolean temporary) { typeIsTemporary = temporary; }
	public boolean cloneAllFields()
	{
		for(int n = 0; n < dbFields.length; n++)
			dbFields[n].setCloneFlag(true);
		return true;
	}

    /**
	 * copyArrayValues
	 *
	 * @param arrayObject Object
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyArrayValues(Object arrayObject) throws IOException, IllegalAccessException
	{
		int size = Array.getLength(arrayObject);
		Object newArray = Array.newInstance(getClassType(), size);
		for(int n = 0; n < size; n++)
		{
			Object oldObject = Array.get(arrayObject, n);
			if(oldObject == null) continue;
			StsDBTypeObject elementClassType = (StsDBTypeObject)StsObjectCopier.getDBTypeList().getDBType(oldObject);
			Object newObject = elementClassType.copyObject(oldObject);
			Array.set(newArray, n, newObject);
		}
		return newArray;
	}

	/**
	 * copyField
	 *
	 * @param oldObject Object
	 * @param newObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 */
	public void copyField(Object oldObject, Object newObject, Field field) throws IOException, IllegalAccessException
	{
		Object oldFieldObject = field.get(oldObject);
		Object newFieldObject = copyObject(oldFieldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(newObject, newFieldObject);
	}

	protected void copyFields(Object oldObject, Object newObject, boolean skipIndexField)
	{
		try
		{
			setFieldsFromTreeMap();
			if (dbFields == null) return;
			int nFields = dbFields.length;
			if (debug) System.out.print("copying " + nFields + " fields for object of class: " + getTypeName() + "  ");
			for (int i = 0; i < nFields; i++)
			{
				StsDBField dbField = (StsDBField) dbFields[i];
				if (debug) System.out.print("  " + i + ": " + dbField.getFieldName());
				if (dbField.getCloneStsObject())
					dbField.clone(oldObject, newObject);
				else
				{
					if(!skipIndexField  || !dbField.getFieldName().equals("index"))
						dbField.copy(oldObject, newObject);
				}
			}
			if(debug) System.out.println();
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "copyObject(in, obj)", e);
		}
	}

	public void cloneField(Object oldObject, Object newObject, Field field) throws IllegalAccessException
	{
		Object fieldObject = field.get(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(newObject, fieldObject);
	}

	/**
	 * readArrayValues
	 *
	 * @param in StsDBInputStream
	 * @return Object
	 * @throws IOException
	 */
	public Object readArrayValues(StsDBInputStream in, Object arrayObject) throws IOException
	{
        Object[] array = (Object[])arrayObject;
        StsDBTypeObject dbType = null;
        int n = 0;
        int size = array.length;
        int typeIndex = -1;
        try
        {
            for(n = 0; n < size; n++)
            {
                typeIndex = in.readInt();
                if (typeIndex == NULL_REFERENCE)
                {
                    array[n] = null;
                    continue;
                }
                if (typeIndex == -1)
                {
                    System.out.println("Should never be the case");
                }
                dbType = (StsDBTypeObject)in.getInputDBType(typeIndex);
                array[n] = dbType.readObject(in);
            }
            return array;
        }
        catch(Exception e)
        {
            String typename = "null";
            if(dbType != null) typename = dbType.getTypeName();
            StsException.outputWarningException(this, "readArrayValues", " for object array[" + n + "]" + " of type: " + typename + " typeIndex: " + typeIndex, e);
            return null;
        }
    }

	/**
	 * readField
	 *
	 * @param in StsDBInputStream
	 * @param classObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
	public void readField(StsDBInputStream in, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		int classListIndex = in.readInt();
		if (classListIndex == NULL_REFERENCE)
		{
			// A null field indicates that the current definition no longer has that field
			if (field != null)
				field.set(classObject, null);
			return;
		}

		StsDBTypeObject dbClassType =  (StsDBTypeObject)in.getInputDBType(classListIndex);
		if (dbClassType == null)
		{
			throw new NullPointerException("StsDBTypeClass::readField - dbClassType cannot be null");
		}
		// The written field may not be the same as the field reference type - could be a sub-class.
		// Therefore, use the actual type to read.
		Object object = null;
		if (field != null)
		{
			object = field.get(classObject);
		}
		if (object == null)
		{
			object =  dbClassType.readObject(in);
		}
		else
		{
			object = dbClassType.readObject(in, object);
		}
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(classObject, object);
	}

	/**
	 * writeArrayValues
	 *
	 * @param out StsDBOutputStream
	 * @param arrayObject Object
	 * @throws IOException
	 */
	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException, IllegalAccessException
	{
		int size = Array.getLength(arrayObject);
		for (int n = 0; n < size; n++)
		{
			Object value = Array.get(arrayObject, n);
			if (value == null)
			{
				out.writeInt(NULL_REFERENCE);
				continue;
			}
			StsDBType dbt = out.getOutputDBType(value);
			StsDBTypeObject dbClassType = (StsDBTypeObject)dbt;
			int classIndex = dbClassType.getIndex();
			if (classIndex == -1)
			{
				System.out.println("Should never get here.");
			}
			out.writeInt(classIndex);
			dbClassType.writeObject(out, value);
		}
	}

	public void exportArrayValues(StsDBFileObjectTrader objectTrader, Object arrayObject) throws IllegalAccessException
	{
		int size = Array.getLength(arrayObject);
		for (int n = 0; n < size; n++)
		{
			Object value = Array.get(arrayObject, n);
			if (value == null)
			{
				continue;
			}
			StsDBType dbt = objectTrader.getCurrentDBType(value);
			StsDBTypeObject dbClassType = (StsDBTypeObject)dbt;
			int classIndex = dbClassType.getIndex();
			if (classIndex == -1)
			{
				System.out.println("Should never get here.");
			}
			dbClassType.exportObject(objectTrader, value);
		}
	}

	/**
	 * writeField
	 *
	 * @param out StsDBOutputStream
	 * @param classObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
    //TODO need if this class is a secondary object (not referenceable),
    //TODO we need to check for circular references to this object which may be in member levels below it
    //TODO i.e., a member object may itself have a member reference to the parent
	public void writeField(StsDBOutputStream out, Object classObject, Field field) throws IllegalAccessException, IOException
	{
        try
        {
            Object value = field.get(classObject);
            if (value == null)
            {
                out.writeInt(NULL_REFERENCE);
                return;
            }

            Class fieldObjectClass = value.getClass();
            Class objectClass = classObject.getClass();
            if(fieldObjectClass == objectClass)
            {
                StsException.systemError("object and field make circular reference to same class: " + StsToolkit.getSimpleClassname(objectClass));
                return;
            }
            out.writeInt(getIndex());
            writeObject(out, value);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "writeField", e);
        }
	}

	public void exportField(StsDBFileObjectTrader objectTrader, Object classObject, Field field) throws IllegalAccessException
	{
		Object value = field.get(classObject);
		if (value == null)
		{
			return;
		}

		exportObject(objectTrader, value);
	}

	protected void initializeAfterLoad(StsDBObjectTypeList typeList) throws ClassNotFoundException
	{
		initializeFieldTypes(typeList);
		setTreeMapFromFields();
		super.initializeAfterLoad(typeList);
	}

	protected void initializeFieldTypes(StsDBObjectTypeList typeList)
	{
		super.initializeFieldTypes(typeList);
		if(dbFields == null) return;
		for (int i = 0; i < dbFields.length; i++)
		{
			StsDBField f = dbFields[i];
			f.initializeFieldTypes(typeList, classType);
		}
	}

	private void setTreeMapFromFields()
	{
		dbFieldsMap = new TreeMap();
		if(dbFields == null) return;
		for (int i = 0; i < dbFields.length; i++)
		{
			StsDBField f = (StsDBField)dbFields[i];
			dbFieldsMap.put(f.getFieldName(), f);
		}
	}

	private void setFieldsFromTreeMap()
	{
		if(dbFields != null) return;

		if(dbFieldsMap == null)
			dbFields = new StsDBField[0];
		else
		{
			Object[] dbFieldObjects = dbFieldsMap.values().toArray();
			int nFields = dbFieldObjects.length;
			dbFields = new StsDBField[nFields];
			for(int n = 0; n < nFields; n++)
				dbFields[n] = (StsDBField)dbFieldObjects[n];
		}
	}

	/* get all the fields for this class including it's superclasses */
	protected TreeMap getAllFields()
	{
		TreeMap sortedFields = new TreeMap();
		try
		{
			// get all class and superclass fields in a TreeMap (key-value pairs sorted by key name)
			for (Class nextClass = getClassType(); nextClass != null; nextClass = nextClass.getSuperclass() )
			{
				StsDBField[] fields = getDBFields(nextClass);
				if (fields == null || fields.length == 0) continue;
				int nFields = fields.length;
				for (int n = 0; n < nFields; n++)
				{
					StsDBField dbField = fields[n];
					String fieldName = dbField.getFieldName();
					StsDBField matchingDBField = (StsDBField)sortedFields.get(fieldName);
					if (matchingDBField != null)
					{
						StsException.systemError(this, "getAllFields", " found duplicate field name " + fieldName +
												 " in classes: " + matchingDBField.getField().getDeclaringClass().getName() +
												 " and " + dbField.getField().getDeclaringClass().getName() + "\n" +
												 "  Will ignore new field found.");
					}
					else
					{
						sortedFields.put(fieldName, dbField);
					}
				}
			}
			return sortedFields;
		}
		catch(Exception e)
		{
			StsException.outputException("DBClass.getFields() failed.", e, StsException.WARNING);
			return null;
		}
	}

	/**
	 *  Gets all declared fields for this class
	 */
	protected StsDBField[] getDBFields(Class c)
	{
		StsDBField[] fields = null;
		int numFields = 0;

		try
		{
			Field[] allFields = c.getDeclaredFields();
			// printFields(allFields, "DBClass.getFields().allFields for class: " + c.getName());
			StsDBField[] tempFields = new StsDBField[allFields.length];
			for (int i = 0; i < allFields.length; i++)
			{
				StsDBField field = StsDBField.constructor(allFields[i]);
				if (field != null)
				{
					tempFields[numFields] = field;
					numFields++;
				}
			}
			// create a list of exact size and copy from the temp list
			fields = new StsDBField[numFields];
			System.arraycopy(tempFields, 0, fields, 0, numFields);
			return fields;
		}
		catch (Exception e)
		{
			StsException.outputException("DBClass.getFields(Class c) failed: ",
										 e, StsException.WARNING);
			return null;
		}
	}

	public boolean representsStsSerialisableObject()
	{
		return true;
	}

	public boolean getMatchesOnDiskVersion()
	{
		return matchesOnDiskVersion;
	}

	public void flagAsWrittenToDisk()
	{
		matchesOnDiskVersion = true;
	}

	protected Object newInstance()
	{
        Constructor constructor;
        Class classType = this.getClassType();
		try
		{
            constructor = classType.getConstructor(new Class[0]);
            Object object = constructor.newInstance();
			if (debug)
			{
				sanityCheck(object);
			}
			return object;
        }
		catch(Exception e)
		{
			StsException.outputWarningException(this, "newInstance", " for " + getTypeName() + ". Probably missing the required null argument constructor for this class.", e);
			return null;
		}
	}

	protected StsDBField[] getDBFields()
	{
		return dbFields;
	}

	public void sanityCheck(Object object)
	{
		StsDBField[] dbFields = getDBFields();
		int nFields = dbFields.length;
		for (int i = 0; i < nFields; i++)
		{
			StsDBField dbField = dbFields[i];
			if (dbField.getDBType().getClass() == StsDBTypeStsClass.class)
			{
				try
				{
					Object fieldObject = dbField.getField().get(object);
					if (fieldObject != null)
					{
						System.out.println("Object sanity check -- Object of type " + object.getClass() + " has been created with a null constructor and");
						System.out.println("   contains one or more non-null references to an StsObject. Check constructor.");
						return;
					}
				}
				catch (IllegalAccessException ex)
				{
				}
				catch (IllegalArgumentException ex)
				{
				}
			}
		}
	}
}
