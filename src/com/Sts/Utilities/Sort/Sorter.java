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

//
//
// Sorter
//
//
public abstract class Sorter 
{
    public void Sort
        (
        Object [] list,
        SortTool tool
        )
    {
        Sort(list,tool,false);
    }

    public abstract void Sort
        (
        Object [] list,
        SortTool tool,
        boolean descending
        );
}

