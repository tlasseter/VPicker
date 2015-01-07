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
// IntegerSortTool
//
//
public class IntegerSortTool 
    implements SortTool
{
    // compare two values
    public int compare
        (
        Object x1,
        Object x2
        )
    {
        if ((x1 instanceof Integer)
        &&  (x2 instanceof Integer))
        {
            int n1 = ((Integer)x1).intValue();
            int n2 = ((Integer)x2).intValue();

            if (n1 < n2)
                return COMP_LESS;
            else
                if (n1 > n2)
                    return COMP_GRTR;
                else
                    return COMP_EQUAL;
        }
        else
            throw SortTool.err1;
    }

}

