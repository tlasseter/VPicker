package com.Sts.DBTypes;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.DB.*;

import java.io.*;

public class StsListItem extends StsSerialize implements Cloneable, StsSerializable, Serializable, Comparable
{
    public int index = -1;
    transient public String name = null;

    public StsListItem()
    {
    }

    public StsListItem(int index, String name)
    {
        this.index = index;
        this.name = name;
    }

    public String toString() { return name; }

    public int compareTo(Object other)
    {
        if(!(other instanceof StsListItem)) return -1;
        StsListItem otherItem = (StsListItem)other;
        return index - otherItem.index;
    }
}
