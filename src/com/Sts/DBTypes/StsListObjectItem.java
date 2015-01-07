package com.Sts.DBTypes;

import com.Sts.DB.*;

import java.io.*;

public class StsListObjectItem extends StsSerialize implements Cloneable, StsSerializable, Serializable, Comparable
{
    public int index = -1;
    transient public Object object;
    transient public String name;

    public StsListObjectItem()
    {
    }

    public StsListObjectItem(int index, String name, Object object)
    {
        this.index = index;
        this.name = name;
        this.object = object;
    }

    public String toString() { return name; }

    public int compareTo(Object other)
    {
        if(!(other instanceof StsListObjectItem)) return -1;
        StsListObjectItem otherItem = (StsListObjectItem)other;
        return index - otherItem.index;
    }
}
