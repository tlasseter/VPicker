package com.Sts.IO;

import java.io.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public interface StsFilenameFilterFace extends FilenameFilter
{

	public String getFilenameEnding(String name);
    public String getFilenameName(String name);
}
