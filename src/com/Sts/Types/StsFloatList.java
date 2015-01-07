
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;

import java.io.*;

public class StsFloatList implements Serializable
{
    private int size = 0;
    private int incSize;
	private int numberDeletedValues = 0;

   	public float[] list;

    static public final int NO_MATCH = StsParameters.NO_MATCH;

	public StsFloatList(int initSize, int incSize)
	{
        this.incSize = incSize;
        initList(initSize);
	}

	public StsFloatList(int initSize)
	{
        incSize = 1;
        initList(initSize);
    }

  	private void initList(int initSize)
    {
    	list = new float[initSize];
    }

    public void add(float value)
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
            StsException.outputException(new StsException("StsFloatList.add(int)",
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
    	checkIndex(index, "StsFloatList.delete(int)");
    	list[index] = -1.0f;
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

    public int getIndex(float value)
    {
    	for (int n = 0; n < list.length; n++)
        	if (list[n] == value) return n;

        return NO_MATCH;
    }

    protected void resize(int newSize)
    {
        float[] oldList = list;
        list = new float[newSize];
        System.arraycopy(oldList, 0, list, 0, size);
   	}

   	public int getSize()
   	{
   		return size;
   	}

    public float getElement(int index)
    {
    	checkIndex(index, "StsFloatList.getElement(int)");
    	return list[index];
    }

    public float getFirst()
    {
    	return getElement(0);
    }

    public float getLast()
    {
    	return getElement(size-1);
    }

    public void deleteLast()
    {
    	delete(size-1);
    }

	public void insertBefore(int index, float newValue)
    {
        if (list.length <= size) resize(list.length+incSize);
        size++;

        for(int n = list.length-1; n > index; n--)
        	list[n] = list[n-1];

        list[index] = newValue;
    }

	public void insertBefore(float oldValue, float newValue)
    {
        int index = getIndex(oldValue);

        if(index == NO_MATCH)
        {
            StsException.outputException( new StsException(StsException.WARNING,
                "StsList.insertBefore(Object, newObject) failed.",
                "Couldn't find matching Object in list."));
        }
        insertBefore(index, newValue);
    }

    public void insertAfter(int index, float newValue)
    {
        if (list.length <= size) resize(list.length+incSize);
        size++;

        for(int n = list.length-1; n > index+1; n--)
        	list[n] = list[n-1];

        list[index+1] = newValue;
    }

    public void unInsert(float value)
    {
        int index = getIndex(value);
        for(int i=index; i<list.length-1; i++)
        	list[i] = list[i+1];
        size--;
    }

	public void insertAfter(float oldValue, float newValue)
    {
        int index = getIndex(oldValue);

        if(index == NO_MATCH)
        {
            StsException.outputException( new StsException(StsException.WARNING,
                "StsList.StsList.insertAfter(oldFloatValue, newFloatValue) failed.",
                "Couldn't find value in list matching: " + oldValue));
        }
        insertAfter(index, newValue);
    }

 	public void addList(StsFloatList extraList)
    {
    	int extraSize = extraList.getSize();
       	if(list.length < size + extraSize) resize(size+extraSize+incSize);

        for(int n = 0; n < extraSize; n++)
         	add(extraList.list[n]);
    }

    public void debugPrint()
    {
    	System.out.println("StsFloatList mainDebug print.  list.length=" + list.length +
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
    	StsFloatList intList = null;

        intList = new StsFloatList(10, 10);

       	int i = 0;

        for(int n = 0; n < 20; n++)
        {
            intList.add((float)i);
            i++;
        }

        intList.debugPrint();
	}
}
