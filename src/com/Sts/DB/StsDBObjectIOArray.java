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

public class StsDBObjectIOArray implements StsDBObjectIO
{
	private StsDBTypeArray arrayType;

	public StsDBObjectIOArray(StsDBTypeObject classType)
	{
		this.arrayType = (StsDBTypeArray)classType;
	}

	/**
	 * copyObject
	 *
	 * @param oldObject Object
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyObject(Object object) throws IOException, IllegalAccessException
	{
		StsDBObjectTypeList typeList = StsObjectCopier.getDBTypeList();
		//JKF don't like the look of this....
		arrayType.setClassFromObject(object);
		arrayType.initialize(typeList);
		return copyArray(object, arrayType.getClassType(), arrayType.getDimension());
	}

	private Object copyArray(Object arrayObject, Class inArrayType, int arrayDimension) throws IOException, IllegalAccessException
	{
		int size = Array.getLength(arrayObject);

		if (arrayDimension == 1)
		{
			return arrayType.copyArrayValues(arrayObject);
		}
		else
		{
			inArrayType = inArrayType.getComponentType();
			Object array = Array.newInstance(inArrayType, size);
			for (int n = 0; n < size; n++)
			{
				Object value = copyArray(array, inArrayType, arrayDimension - 1);
				Array.set(array, n, value);
			}
			return array;
		}
	}

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in) throws IOException
	{
		return readArray(in, arrayType.getClassType(), arrayType.getDimension());
	}

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @param obj Object
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in, Object obj) throws IOException
	{
		return readObject(in);
	}

	private Object readArray(StsDBInputStream in, Class inArrayType, int arrayDimension) throws IOException
	{
		boolean isNull = in.readBoolean();
		if (isNull) return null;

		int size = in.readInt();

		if (arrayDimension == 1)
		{
			Object newArray = Array.newInstance(inArrayType.getComponentType(), size);
			return arrayType.readArrayValues(in, newArray);
		}
		else
		{
			inArrayType = inArrayType.getComponentType();
			Object array = Array.newInstance(inArrayType, size);
			for (int n = 0; n < size; n++)
			{
				Object value = readArray(in, inArrayType, arrayDimension - 1);
				Array.set(array, n, value);
			}
			return array;
		}
	}

	/**
	 * writeObject
	 *
	 * @param out StsDBOutputStream
	 * @param obj Object
	 * @throws IOException
	 * @throws IllegalAccessException
	 */
	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		writeArray(out, obj, arrayType.getDimension());
	}

	private void writeArray(StsDBOutputStream out, Object arrayObject, int arrayDimension) throws IOException, IllegalAccessException
	{
		out.writeBoolean(arrayObject == null);
		if (arrayObject == null) return;

		int size = Array.getLength(arrayObject);
		out.writeInt(size);

		if (arrayDimension == 1)
		{
			arrayType.writeArrayValues(out, arrayObject);
		}
		else
		{
			for (int n = 0; n < size; n++)
			{
				Object subArray = Array.get(arrayObject, n);
				writeArray(out, subArray, arrayDimension - 1);
			}
		}
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
		exportArray(objectTrader, obj, arrayType.getDimension());
	}

	private void exportArray(StsDBFileObjectTrader objectTrader, Object arrayObject, int arrayDimension) throws IllegalAccessException
	{
		if (arrayObject == null) return;

		int size = Array.getLength(arrayObject);

		if (arrayDimension == 1)
		{
			arrayType.exportArrayValues(objectTrader, arrayObject);
		}
		else
		{
			for (int n = 0; n < size; n++)
			{
				Object subArray = Array.get(arrayObject, n);
				exportArray(objectTrader, subArray, arrayDimension - 1);
			}
		}
	}
}
