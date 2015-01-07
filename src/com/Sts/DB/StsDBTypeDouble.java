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

public class StsDBTypeDouble extends StsDBType implements StsSerializable
{
	public StsDBTypeDouble()
	{
		super(Double.TYPE);
	}

	public StsDBTypeDouble(int index)
	{
		super(index, Double.TYPE);
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
		return (double[])arrayObject;
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
		double value = field.getDouble(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setDouble(newObject, value);
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
		double[] doubleArray = (double[])arrayObject;
		int size = doubleArray.length;
		for(int n = 0; n < size; n++)
		{
			doubleArray[n] = in.readDouble();
		}
		return doubleArray;
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
		double value = in.readDouble();
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setDouble(classObject, value);
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
		double[] doubleArray = (double[])arrayObject;
		int size = doubleArray.length;
		for(int n = 0; n < size; n++)
		{
			out.writeDouble(doubleArray[n]);
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
		double value = field.getDouble(classObject);
		out.writeDouble(value);
	}
}
