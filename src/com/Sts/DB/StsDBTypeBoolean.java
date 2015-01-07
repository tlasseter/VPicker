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

public class StsDBTypeBoolean extends StsDBType implements StsSerializable
{
	public StsDBTypeBoolean()
	{
		super(Boolean.TYPE);
	}

	public StsDBTypeBoolean(int index)
	{
		super(index, Boolean.TYPE);
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
		return (boolean[])arrayObject;
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
		boolean value = field.getBoolean(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setBoolean(newObject, value);
	}

	/**
	 * readArrayValues
	 *
	 * @param in StsDBInputStream
	 * @param size int
	 * @return Object
	 * @throws IOException
	 */
	public Object readArrayValues(StsDBInputStream in, Object arrayObject) throws IOException
	{
		boolean[] booleanArray = (boolean[])arrayObject;
		int size = booleanArray.length;
		for(int n = 0; n < size; n++)
		{
			booleanArray[n] = in.readBoolean();
		}
		return booleanArray;
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
		boolean value = in.readBoolean();
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setBoolean(classObject, value);
	}

	/**
	 * writeArrayValues
	 *
	 * @param out StsDBOutputStream
	 * @param arrayObject Object
	 * @param size int
	 * @throws IOException
	 */
	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException
	{
		boolean[] booleanArray = (boolean[])arrayObject;
		int size = booleanArray.length;
		for(int n = 0; n < size; n++)
		{
			out.writeBoolean(booleanArray[n]);
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
		boolean value = field.getBoolean(classObject);
		out.writeBoolean(value);
	}
}
