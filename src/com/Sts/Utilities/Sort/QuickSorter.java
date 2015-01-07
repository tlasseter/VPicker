//  JAVA ALGORITHMS
//  ---------------
//  Copyright 1997 Scott Robert Ladd
//  All rights reserved
//
//  This software source code is sold as a component of the
//  book JAVA ALGORITHMS, written by Scott Robert Ladd and
//  published by McGraw-Hill, Inc. Please read the LICENSE
//  AGREEMENT and DISCLAIMER OF WARRANTY printed in the book.
//
//  You may freely compile and link this source code into your 
//  non-commercial software programs, providing that you do
//  not redistribute the source code or object code derived
//  therefrom. If you want to use this source code in a
//  commercial application, you must obtain written permission
//  by contacting:
//
//      Scott Robert Ladd
//      P.O. Box 617
//      Silverton, Colorado
//      81433-0617 USA
//
//  This software is sold "as is" without warranty of any kind.

package com.Sts.Utilities.Sort;

import com.Sts.Types.*;
//
//
// QuickSorter
//
//
public class QuickSorter 
    extends Sorter
{
    private IllegalArgumentException err1 =
        new IllegalArgumentException("stack overflow in QuickSort");

    private class StackItem
    {
        public int left;
        public int right;
    }

    public void Sort
        (
        Object [] list,
        SortTool tool,
        boolean descending
        )
    {
        // create stack
        final int stackSize = 32;
        StackItem [] stack = new StackItem [stackSize];

        for (int n = 0; n < 32; ++n)
            stack[n] = new StackItem();

        int stackPtr = 0;

        // determine direction of sort
        int comp;

        if (descending)
            comp = SortTool.COMP_GRTR;
        else
            comp = SortTool.COMP_LESS;

        // size of minimum partition to median-of-three
        final int Threshold = 7;

        // sizes of left and right partitions
        int lsize, rsize;

        // create working indexes
        int l, r, mid, scanl, scanr, pivot;

        // get the list length by finding the number of non-nulls
        int length = 0;
        while(length < list.length && list[length] != null) length++;

        l = 0;
        r = length - 1;

        Object temp;
    
        // main loop
        while (true)
        {
            while (r > l)
            {
                if ((r - l) > Threshold)
                {
                    // "median-of-three" partitioning
                    mid = (l + r) / 2;
    
                    // three-sort left, middle, and right elements
                    if (tool.compare(list[mid],list[l]) == comp)
                    {
                        temp      = list[mid];
                        list[mid] = list[l];
                        list[l]   = temp;
                    }

                    if (tool.compare(list[r],list[l]) == comp)
                    {
                        temp    = list[r];
                        list[r] = list[l];
                        list[l] = temp;
                    }

                    // three-sort left, middle, and right elements
                    if (tool.compare(list[r],list[mid]) == comp)
                    {
                        temp      = list[mid];
                        list[mid] = list[r];
                        list[r]   = temp;
                    }

                    // set-up for partitioning
                    pivot = r - 1;
    
                    temp        = list[mid];
                    list[mid]   = list[pivot];
                    list[pivot] = temp;
    
                    scanl = l + 1;
                    scanr = r - 2;
                }
                else
                {
                    // set-up for partitioning
                    pivot = r;
                    scanl = l;
                    scanr = r - 1;
                }
    
                for (;;)
                {
                    // scan from left for element >= to pivot
                    while ((tool.compare(list[scanl],list[pivot]) == comp) && (scanl < r))
                        ++scanl;
    
                    // scan from right for element <= to pivot
                    while ((tool.compare(list[pivot],list[scanr]) == comp) && (scanr > l))
                        --scanr;
    
                    // if scans have met, exit inner loop
                    if (scanl >= scanr)
                        break;
    
                    // exchange elements
                    temp        = list[scanl];
                    list[scanl] = list[scanr];
                    list[scanr] = temp;
    
                    if (scanl < r)
                        ++scanl;
    
                    if (scanr > l)
                        --scanr;
                }
    
                // exchange final element
                temp        = list[scanl];
                list[scanl] = list[pivot];
                list[pivot] = temp;
    
                // place largest partition on stack
                lsize = scanl - l;
                rsize = r - scanl;
    
                if (lsize > rsize)
                {
                    if (lsize != 1)
                    {
                        ++stackPtr;
    
                        if (stackPtr == stackSize)
                            throw err1;
    
                        stack[stackPtr].left  = l;
                        stack[stackPtr].right = scanl - 1;
                    }
    
                    if (rsize != 0)
                        l = scanl + 1;
                    else
                        break;
                }
                else
                {
                    if (rsize != 1)
                    {
                        ++stackPtr;
    
                        if (stackPtr == stackSize)
					    	throw err1;
    
                        stack[stackPtr].left  = scanl + 1;
                        stack[stackPtr].right = r;
                    }
    
                    if (lsize != 0)
                        r = scanl - 1;
                    else
                        break;
                }
            }
    
            // iterate with values from stack
            if (stackPtr != 0)
            {
                l = stack[stackPtr].left;
                r = stack[stackPtr].right;
    
                --stackPtr;
            }
            else
                break;
        }
    }

    public static void main(String[] args)
    {
        int i, ii;

        StsList list = new StsList(5, 2);

        System.out.println("Input values: ");
        for (i = 0; i < 10; ++i)
        {
            ii = (int)(100*Math.random());
            System.out.print(" " + ii);
            list.add(new Integer(ii));
        }

        Sorter sort = new QuickSorter();
        IntegerSortTool tool = new IntegerSortTool();
        sort.Sort(list.list, tool, false);

        System.out.println("Output values: ");
        int size = list.getSize();
        for(i = 0; i < size; i++)
        {
            Integer intObj = (Integer)list.list[i];
            ii = intObj.intValue();
//            ii = (Integer)(list.list[i]).intValue();
            System.out.print(" " + ii);
        }
	}
}

