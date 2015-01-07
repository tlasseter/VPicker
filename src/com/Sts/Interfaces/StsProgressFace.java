package com.Sts.Interfaces;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public interface StsProgressFace
{
    public static final int WARNING = 0;
    public static final int INFO = 1;
    public static final int FLASH = 2;
    public static final int FATAL = 3;
    public static final int ERROR = 4;

    public void appendLine(String line);
    public void initialize(int maxValue);
    public void setMaximum(int maxValue);
    public void setValue(int value);
    public void setDescriptionAndLevel(String progressDescription, int level);
    public void cancel();
}
