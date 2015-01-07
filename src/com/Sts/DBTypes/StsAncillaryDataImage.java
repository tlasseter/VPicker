package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsAncillaryDataImage extends StsAncillaryData implements StsTreeObjectI
{
    static StsObjectPanel objectPanel = null;

    public static String dataType = "Image";

    public static byte CGM = 0;
    public static byte JPG = 1;
    public static byte BMP = 2;
    public static byte GIF = 3;
    public static byte TIF = 4;

    public static String[] EXTENSIONS = {"cgm", "jpg", "bmp", "gif", "tif"};
    public static byte types[] = { CGM, JPG, BMP, GIF, TIF };
    public static String[] msCommands = { "", "mspaint.exe", "mspaint.exe", "mspaint.exe", "mspaint.exe" };
    public static String[] unixCommands = { "", "", "", "", "" };

    public StsAncillaryDataImage()
    {
    }

    public StsAncillaryDataImage(String filename, String cmd, StsColor color, boolean persist)
    {
        super(filename, cmd, color, currentModel.getProject().getXOrigin(), currentModel.getProject().getYOrigin(),
             (double)currentModel.getProject().getZorTMin(), persist);
    }

    public StsAncillaryDataImage(String filename, StsColor color, boolean persist)
    {
        super(filename, getDefaultCommandFromFilename(filename), color, currentModel.getProject().getXOrigin(), currentModel.getProject().getYOrigin(),
             (double)currentModel.getProject().getZorTMin(), persist);
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsFieldBean[] getDisplayFields()
    {
        return displayFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }
    static public String getFileExtString(byte type) { return EXTENSIONS[type]; }
    static public String[] getFileExtensions() { return EXTENSIONS; }
    static public String getDefaultCommandFromFilename(String filename)
    {
        for(int i=0; i<EXTENSIONS.length; i++)
        {
            if(filename.toLowerCase().endsWith(EXTENSIONS[i]))
            {
                String osName = System.getProperty("os.name").toString();
                if(osName.indexOf("Windows") >= 0)
                    return msCommands[i];
                else
                    return unixCommands[i];
            }
        }
        return null;
    }
    static public boolean isThisType(String filename)
    {
        for(int i=0; i<EXTENSIONS.length; i++)
        {
            if(filename.toLowerCase().endsWith(EXTENSIONS[i]))
            {
                return true;
            }
        }
        return false;
    }
}
