package com.Sts.DB;

import java.io.*;

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

public interface StsCustomSerializable extends StsSerializable
{
	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException;
	public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException;
    public void exportObject(StsDBFileObjectTrader objectTrader) throws IllegalAccessException;
}
