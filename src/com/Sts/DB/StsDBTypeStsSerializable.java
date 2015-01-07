package com.Sts.DB;

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

public class StsDBTypeStsSerializable extends StsDBTypeClass implements StsSerializable
{
	public StsDBTypeStsSerializable()
	{
		super();
	}

	public StsDBTypeStsSerializable(int index, Class classType)
	{
		super(index, classType);
	}

	/* get all the fields for this class including it's superclasses */
	protected TreeMap getAllFields()
	{
		TreeMap sortedFields = new TreeMap();
		try
		{
			// get all class and superclass fields in a TreeMap (key-value pairs sorted by key name)
			Class classType = getClassType();
			for (Class nextClass = classType; nextClass != null; nextClass = nextClass.getSuperclass())
			{
				if (! StsSerializable.class.isAssignableFrom(nextClass))
				{
					return sortedFields;
				}
				StsDBField[] fields = getDBFields(nextClass);
				if( fields == null || fields.length == 0 ) continue;
				int nFields = fields.length;
				for(int n = 0; n < nFields; n++)
				{
					StsDBField dbField = fields[n];
					String fieldName = dbField.getFieldName();
					StsDBField matchingDBField = (StsDBField)sortedFields.get(fieldName);
					if (matchingDBField != null)
					{
						StsException.systemError(this, "getAllFields",  " found duplicate field name " + fieldName +
														 " for classtype " + classType.getName() + " in superClass " + nextClass.getName() + "\n" +
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
}
