
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsIntList implements Serializable
{
    private int size = 0;
    private int incSize;
	private int numberDeletedValues = 0;

   	public int[] list;

	public StsIntList(int initSize, int incSize)
	{
        this.incSize = incSize;
        initList(initSize);
	}

	public StsIntList(int initSize)
	{
        incSize = 1;
        initList(initSize);
    }

  	private void initList(int initSize)
    {
    	list = new int[initSize];
    }

    public void add(int value)
    {
        if(numberDeletedValues > 0)
        {
        	for(int n = 0; n < size; n++)
            {
            	if(list[n] == -1)
                {
                	list[n] = value;
                    numberDeletedValues--;
                    return;
                }
            }
            StsException.outputException(new StsException("StsIntList.add(int)",
                "System error: didn't find deleted entry in list."));
        }
        else
        {
       		if(list.length == size) resize(size+incSize);
        	list[size++] = value;
        }
    }

    public void delete(int index)
    {
    	checkIndex(index, "StsIntList.delete(int)");
    	list[index] = -1;
        numberDeletedValues++;
    }

    private void checkIndex(int index, String methodDescription)
    {
    	if(index+1 > size)
            StsException.outputException(new StsException(methodDescription,
                "Specified index: ", index, " > index of last item: ", size-1));
        if(index < 0)
            StsException.outputException(new StsException(methodDescription,
                "Specified index: ", index, " <  ", 0));
    }

    public void add(StsObject obj)
    {
        if(list.length == size) resize(size+incSize);
        list[size++] = obj.getIndex();
    }

    public void add(StsObject obj, int sign)
    {
        if(list.length == size) resize(size+incSize);

        if(sign >= 0)
        	list[size++] = obj.getIndex();
        else
        	list[size++] = -obj.getIndex();
    }

    protected void resize(int newSize)
    {
        int[] oldList = list;
        list = new int[newSize];
        System.arraycopy(oldList, 0, list, 0, size);
   	}

   	public int getSize()
   	{
   		return size;
   	}

    public int getElement(int index)
    {
    	checkIndex(index, "StsIntList.getElement(int)");
    	return list[index];
    }

    public int getFirst()
    {
    	return getElement(0);
    }

    public int getLast()
    {
    	return getElement(size-1);
    }

    public void deleteLast()
    {
    	delete(size-1);
    }


 	public void addList(StsIntList extraList)
    {
    	int extraSize = extraList.getSize();
       	if(list.length < size + extraSize) resize(size+extraSize+incSize);

        for(int n = 0; n < extraSize; n++)
         	add(extraList.list[n]);
    }

    public void debugPrint()
    {
    	System.out.println("StsIntList mainDebug print.  list.length=" + list.length +
        				   " intList.size=" + size);

    	if(size > 20)
        {
    		System.out.print("First 10 entries: ");
        	for(int n = 0; n < 10; n++)
        		System.out.print(" " + list[n]);

            System.out.println("");

            System.out.print("Last 10 entries: ");
            for(int n = size-10; n < size; n++)
        		System.out.print(" " + list[n]);
		}
        else
        {
            System.out.print("List entries: ");

        	for(int n = 0; n < size; n++)
            	System.out.print(" " + list[n]);
        }
    }

	public static void main(String[] args)
	{
    	StsIntList intList = null;

        intList = new StsIntList(10, 10);

       	int i = 0;

        for(int n = 0; n < 20; n++)
        {
            intList.add(i);
            i++;
        }

        intList.debugPrint();
	}
}

