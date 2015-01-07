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

public class StsDBTypeLong extends StsDBType implements StsSerializable
{
	public StsDBTypeLong()
	{
		super(Long.TYPE);
	}

	public StsDBTypeLong(int index)
	{
		super(index, Long.TYPE);
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
		return (long[])arrayObject;
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
		long value = field.getLong(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setLong(newObject, value);
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
		long[] longArray = (long[])arrayObject;
		int size = longArray.length;
		for(int n = 0; n < size; n++)
		{
			longArray[n] = in.readLong();
		}
		return longArray;
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
		long value = in.readLong();
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setLong(classObject, value);
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
		long[] longArray = (long[])arrayObject;
		int size = longArray.length;
		for(int n = 0; n < size; n++)
		{
			out.writeLong(longArray[n]);
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
		long value = field.getLong(classObject);
		out.writeLong(value);
	}
}
