package com.Sts.DB;

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

public class StsDBTypePrimitiveObject extends StsDBTypeObject implements StsSerializable
{
	public StsDBTypePrimitiveObject()
	{
		super();
	}

	public StsDBTypePrimitiveObject(int index, Class classType, StsDBObjectIO objectIO)
	{
		super(index, classType, objectIO);
	}

	protected StsDBObjectIO getObjectIO()
	{
		throw new RuntimeException("StsDBPrimitiveObjectType::initiliseObjectIO - should never get here!");
	}

	/**
	 * copyArrayValues
	 *
	 * @param arrayObject Object
	 * @param size int
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyArrayValues(Object arrayObject) throws IllegalAccessException
	{
		return arrayObject;
	}

	/**
	 * copyField
	 *
	 * @param oldObject Object
	 * @param newObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 */
	public void copyField(Object oldObject, Object newObject, Field field) throws IllegalAccessException
	{
		Object object = field.get(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(newObject, object);
	}

	/**
	 * readArrayValues
	 *
	 * @param in StsDBInputStream
	 * @param arrayObject Object
	 * @return Object
	 * @throws IOException
	 */
	public Object readArrayValues(StsDBInputStream in, Object arrayObject) throws IOException
	{
		Object[] objectArray = (Object[])arrayObject;
		int size = objectArray.length;
		for(int n = 0; n < size; n++)
		{
			objectArray[n] = readObject(in);
		}
		return objectArray;
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
		Object object = readObject(in);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(classObject, object);
	}

	/**
	 * writeArrayValues
	 *
	 * @param out StsDBOutputStream
	 * @param arrayObject Object
	 * @param size int
	 * @throws IOException
	 * @throws IllegalAccessException
	 */
	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException, IllegalAccessException
	{
		Object[] objectArray = (Object[])arrayObject;
		int size = objectArray.length;
		for(int n = 0; n < size; n++)
		{
			writeObject(out, objectArray[n]);
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
	public void writeField(StsDBOutputStream out, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		Object fieldObject = field.get(classObject);
		// If the field type is an object or Number and a primitive object is assigned then we need to add additional
		// type information so that it can be read in correctly
		Class fieldType = field.getType();
		if (fieldObject != null && (fieldType == Object.class || fieldType == java.lang.Number.class))
		{
			out.writeInt(getIndex());
		}
		writeObject(out, fieldObject);
	}

	protected Object newInstance()
	{
		throw new RuntimeException("Cannot instantiate a primitive object");
	}
}
