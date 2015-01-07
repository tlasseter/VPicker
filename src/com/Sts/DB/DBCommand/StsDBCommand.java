package com.Sts.DB.DBCommand;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;

public abstract class StsDBCommand
{
	static final protected transient boolean debug = Main.isDbCmdDebug;

	// When adding a new class the following steps need to happen....
	// 1. Create a new 'static final byte' index for the command below. This index must be unique.
	// 2. Increment the static int 'NUMBER_OF_COMMANDS'.
	// 3. Add the new command class to the 'dbCommandClasses' array using the new index.
	// 4. Implement the new command and return the new index from the method 'getDBCommandClassIndex'.

	private final static int NUMBER_OF_COMMANDS = 20;

	// These commands indexes should never be modified to maintain backwards compatibility.
	static public final byte ADD_REF_CMD_INDEX = 0;
	static public final byte CHANGE_CMD_INDEX = 1;
	static public final byte DB_METHOD_CMD_INDEX = 2;
	static public final byte DB_TYPE_CMD_INDEX = 3;
	static public final byte DELETE_REF_CMD_INDEX = 4;
	static public final byte EXPORT_ROOT_OBJECT_CMD_INDEX = 5;
	static public final byte INSTANCE_ADD_CMD_INDEX = 6;
	static public final byte INSTANCE_CHG_CMD_INDEX = 7;
	static public final byte INSTANCE_DELETE_CMD_INDEX = 8;
	static public final byte INSTANCE_EXPORT_CMD_INDEX = 9;
	static public final byte SAVE_MODEL_CMD_INDEX = 10;
	static public final byte SIMPLE_OBJECT_ADD_CMD_INDEX = 11;
	static public final byte START_TRANSACTION_CMD_INDEX = 12;
	static public final byte MODEL_DB_REF_CMD_INDEX = 13;
	static public final byte ARRAY_CHANGE_CMD_INDEX = 14;
	static public final byte ARRAY_INSERT_CMD_INDEX = 15;
    static public final byte ARRAY_DELETE_CMD_INDEX = 16;
    static public final byte ARRAY_ELEMENT_CHANGE_CMD_INDEX = 17;
    static public final byte END_DB_CMD_INDEX = 18;
    static public final byte ADD_TRANS_MODEL_OBJECT_CMD_INDEX = 19;
   // Array is required to allow look up of class given a command index.
	static final Class[] dbCommandClasses = new Class[NUMBER_OF_COMMANDS];

	static
	{
		dbCommandClasses[ADD_REF_CMD_INDEX] = StsAddRefCmd.class;
		dbCommandClasses[CHANGE_CMD_INDEX] = StsChangeCmd.class;
		dbCommandClasses[DB_METHOD_CMD_INDEX] = StsDBMethodCmd.class;
		dbCommandClasses[DB_TYPE_CMD_INDEX] = StsDBTypeCmd.class;
		dbCommandClasses[DELETE_REF_CMD_INDEX] = StsDeleteRefCmd.class;
		dbCommandClasses[EXPORT_ROOT_OBJECT_CMD_INDEX] = StsExportRootObjectCmd.class;
		dbCommandClasses[INSTANCE_ADD_CMD_INDEX] = StsInstanceAddCmd.class;
		dbCommandClasses[INSTANCE_CHG_CMD_INDEX] = StsInstanceChgCmd.class;
		dbCommandClasses[INSTANCE_DELETE_CMD_INDEX] = StsInstanceDeleteCmd.class;
		dbCommandClasses[INSTANCE_EXPORT_CMD_INDEX] = StsInstanceExportCmd.class;
		dbCommandClasses[SAVE_MODEL_CMD_INDEX] = StsSaveModelCmd.class;
		dbCommandClasses[SIMPLE_OBJECT_ADD_CMD_INDEX] = StsSimpleObjectAddCmd.class;
		dbCommandClasses[START_TRANSACTION_CMD_INDEX] = StsStartTransactionCmd.class;
//		dbCommandClasses[MODEL_DB_REF_CMD_INDEX] = StsModelDBRefCmd.class;
		dbCommandClasses[ARRAY_CHANGE_CMD_INDEX] = StsArrayChangeCmd.class;
		dbCommandClasses[ARRAY_INSERT_CMD_INDEX] = StsArrayInsertCmd.class;
        dbCommandClasses[ARRAY_DELETE_CMD_INDEX] = StsArrayDeleteCmd.class;
        dbCommandClasses[ARRAY_ELEMENT_CHANGE_CMD_INDEX] = StsArrayElementChangeCmd.class;
        dbCommandClasses[END_DB_CMD_INDEX] = StsEndDBCmd.class;
        dbCommandClasses[ADD_TRANS_MODEL_OBJECT_CMD_INDEX] = StsAddTransientModelObjectCmd.class;
    }

	public StsDBCommand()
	{
	}

	static final public Class getDBCommandClass(byte index)
	{
        if(index < 0 || index > NUMBER_OF_COMMANDS-1)
        {
            StsException.systemError("StsDBCommand.getDBCommandClass() failed. No command for index " + index);
            return null;
        }
        return dbCommandClasses[index];
	}

	public void abort() throws StsException
	{}

	public abstract void write(StsDBOutputStream dbOutputStream) throws IOException;

	public abstract void read(StsDBInputStream dbInputStream) throws IOException;

	public abstract String toDebugString();

	public abstract byte getDBCommandClassIndex();

    public void debugMessageWrite()
    {
        // long position = StsModel.getCurrentModel().getDatabase().getPosition();
        System.out.println(StsToolkit.getSimpleClassname(this) + ".write() " + toDebugString());
    }
    
    public void debugMessageRead()
	{
        // long position = StsModel.getCurrentModel().getDatabase().getPosition();
        System.out.println(StsToolkit.getSimpleClassname(this) + ".read() " + toDebugString());
	}

	static public String toDebugString(StsObject stsObject)
	{
		if (stsObject != null)
			return StsToolkit.getSimpleClassname(stsObject) + "[" + stsObject.getIndex() + "]";
		else
			return "null StsObject";
	}

	static public String toDebugString(Class c, int objectIndex)
	{
		return "StsObject " + StsToolkit.getSimpleClassname(c) + "[" + objectIndex + "]";
	}

	static public String toDebugString(Object object)
	{
		if (object != null)
			return StsToolkit.getSimpleClassname(object);
		else
			return "null object";
	}

	public void writeCmdClassIndex(StsDBOutputStream dbOutputStream) throws IOException
	{
        int index = getDBCommandClassIndex();
        dbOutputStream.writeByte(index);
    }

    public boolean checkFieldName(StsObject obj, String fieldName)
    {
        if(fieldName == null) return false;
        Field field = obj.getField(fieldName);
        if(field != null) return true;
        StsException.systemError(StsChangeCmd.class, "constructor", "Object " + obj.getClass() + " has no field: " + fieldName);
        return false;
    }


    public boolean checkFieldName(Class c, String fieldName)
    {
        Field field = null;

        if(fieldName == null) return false;
        try
        {
            field = c.getField(fieldName);
        }
        catch(Exception e)
        {
            StsException.systemError(this, "checkFieldName", "failed to find field for class " + c.getName() + " field " + fieldName);
            return false;
        }
        return field != null;
    }
}