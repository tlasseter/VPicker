
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;


public class StsCategoricalFacies
{
    // constants
    static public final int INVALID_FACIES = -1;
    static private final int MAX_FACIES = 9;
    static private final String[] FACIES_NAMES = { "A", "B", "C", "D", "E", "F", "G", "H", "I" };

    static public final int NO_MATCH = StsParameters.NO_MATCH;

    // fields
    private int nFacies = 0;
    private StsList[] faciesLists = null;          // lists of lithologies by facies
    private StsList lithologyList = null;          // list of all lithologies
    private Hashtable lithologyFaciesTable = null; // table of lithology-facies pairs

    /** constructTraceAnalyzer */
    public StsCategoricalFacies(int nFacies)
        throws StsException
    {
        if (nFacies < 2) throw new StsException("Too few facies specified.");
        if (nFacies > MAX_FACIES) throw new StsException("Too many facies specified.");
        this.nFacies = nFacies;

        // build facies lists & lithologies hash table
        faciesLists = new StsList[nFacies];
        for (int i=0; i<nFacies; i++) faciesLists[i] = new StsList(1, 1);
        lithologyFaciesTable = new Hashtable(nFacies*2);
        lithologyList = new StsList(nFacies*2, 1);
    }

    /** return array of facies names used in this object instance */
    public String[] getFaciesNames()
    {
        if (nFacies < 2) return null;
        String[] names = new String[nFacies];
        for (int i=0; i<nFacies; i++) names[i] = FACIES_NAMES[i];
        return names;
    }

    /** return array of lithology numbers for a facies */
    public int[] getLithologies(int facies)
    {
        // check to see if the facies is valid
        if (facies < 0 || facies > nFacies-1) return null;

        int nLithologies = faciesLists[facies].getSize();
        if (nLithologies < 1) return null;
        int[] lithologies = new int[nLithologies];
        for (int i=0; i<nLithologies; i++)
        {
            Integer lithologyObject = (Integer)faciesLists[facies].getElement(i);
            lithologies[i] = lithologyObject.intValue();
        }
        Arrays.sort(lithologies);
        return lithologies;
    }

    /** return array of lithology numbers */
    public int[] getAllLithologies()
    {
        int nLithologies = lithologyList.getSize();
        if (nLithologies < 1) return null;
        int[] lithologies = new int[nLithologies];
        for (int i=0; i<nLithologies; i++)
        {
            Integer lithologyObject = (Integer)lithologyList.getElement(i);
            lithologies[i] = lithologyObject.intValue();
        }
        Arrays.sort(lithologies);
        return lithologies;
    }

    /** return array of unassigned lithology numbers */
    public int[] getUnassignedLithologies()
    {
        int nLithologies = lithologyList.getSize();
        if (nLithologies < 1) return null;
        int[] lithologies = new int[nLithologies];
        int nUnassigned = 0;
        for (int i=0; i<nLithologies; i++)
        {
            Integer lithologyObject = (Integer)lithologyList.getElement(i);
            if (lithologyFaciesTable.get(lithologyObject) == null)
            {
                lithologies[nUnassigned] = lithologyObject.intValue();
                nUnassigned++;
            }
        }
        if (nUnassigned == 0) return null;
        int[] unassignedLithologies = new int[nUnassigned];
        System.arraycopy(lithologies, 0, unassignedLithologies, 0, nUnassigned);
        Arrays.sort(unassignedLithologies);
        return unassignedLithologies;
    }

    /** add lithology values to the lithology list */
    public void addLithologies(float[] lithologies)
    {
        if (lithologies == null) return;
        Arrays.sort(lithologies);
        Integer lithologyObject = new Integer((int)lithologies[0]);
        if (!lithologyList.contains(lithologyObject)) lithologyList.add(lithologyObject);
        int last = (int)lithologies[0];
        for (int i=1; i<lithologies.length; i++)
        {
            if (last == lithologies[i]) continue;
            lithologyObject = new Integer((int)lithologies[i]);
            if (!lithologyList.contains(lithologyObject)) lithologyList.add(lithologyObject);
            last = (int)lithologies[i];
        }
    }

    /** get the facies name for a facies */
    static public String getFaciesName(int facies)
    {
        // check to see if the facies is valid
        if (facies < 0 || facies > MAX_FACIES-1) return null;
        return FACIES_NAMES[facies];
    }

    /** get the facies for a facies name */
    static public int getFaciesCategory(String name)
    {
        if (name == null) return INVALID_FACIES;
        switch (name.charAt(0))
        {
            case 'A':   return 0;
            case 'B':   return 1;
            case 'C':   return 2;
            case 'D':   return 3;
            case 'E':   return 4;
            case 'F':   return 5;
            case 'G':   return 6;
            case 'H':   return 7;
            case 'I':   return 8;
        }
        return INVALID_FACIES;
    }

    /** set a lithology-facies pair (either add or move) */
    public boolean set(int facies, int lithology)
    {
        // check to see if the facies is valid
        if (facies < 0 || facies > nFacies-1) return false;

        // check to see if the lithology is valid
        Integer lithologyObject = new Integer(lithology);
        if (!lithologyList.contains(lithologyObject)) return false;

        // see if it's already there
        if (faciesLists[facies].contains(lithologyObject)) return true;

        // remove this lithology if we already have it
        removeLithology(lithology);

        // add the lithology object to the facies list & put to hash table
        faciesLists[facies].add(lithologyObject);
        lithologyFaciesTable.put(lithologyObject, FACIES_NAMES[facies]);

        return true;
    }

    /** add/replace all the facies-lithology pairs for a facies category */
    public boolean addReplace(int facies, int[] lithologies)
    {
        // check to see if the facies is valid
        if (facies < 0 || facies > nFacies-1) return false;

        // check to see if the lithologies are valid
        if (lithologies != null)
        {
            for (int i=0; i<lithologies.length; i++)
            {
                Integer lithologyObject = new Integer(lithologies[i]);
                if (!lithologyList.contains(lithologyObject)) return false;
            }
        }

        // remove previous lithologies in the hash table & any facies list they're in
        removeFacies(facies);

        // any lithologies to set?
        if (lithologies == null) return true;  // okay

        // remove any new lithologies already in a facies list
        for (int i=0; i<lithologies.length; i++) removeLithology(lithologies[i]);

        // add the lithology objects to the facies lists & hash table
        for (int i=0; i<lithologies.length; i++) set(facies, lithologies[i]);

        return true;
    }

    /** remove a facies-lithology pair */
    public boolean removeLithology(int lithology)
    {
        Integer lithologyObject = new Integer(lithology);
        String faciesName = (String)lithologyFaciesTable.get(lithologyObject);
        if (faciesName == null) return false;
        int faciesCategory = getFaciesCategory(faciesName);
        return (faciesLists[faciesCategory].delete(lithologyObject) != NO_MATCH);
    }

    /** remove all lithologies for a facies */
    public boolean removeFacies(int facies)
    {
        // check to see if the facies is valid
        if (facies < 0 || facies > nFacies-1) return false;

        int nLithologies = faciesLists[facies].getSize();
        for (int i=0; i<nLithologies; i++)
        {
            // remove the lithology object from the hash table
            Integer lithologyObject = (Integer)faciesLists[facies].getElement(i);
            lithologyFaciesTable.remove(lithologyObject);
        }
        faciesLists[facies].deleteAll();

        return true;
    }

    /** delete a lithology completely */
    public boolean delete(int lithology)
    {
        // remove the lithology object from the facies list
        Integer lithologyObject = new Integer(lithology);
        if (lithologyList.delete(lithologyObject) == NO_MATCH) return false;
        return removeLithology(lithology);
    }

    /** delete an array of lithologies completely */
    public boolean delete(int[] lithologies)
    {
        if (lithologies == null) return false;
        for (int i=0; i<lithologies.length; i++)
        {
            if (!delete(lithologies[i])) return false;
        }
        return true;
    }

    /** convert an array of lithologies into an array of facies */
    public float[] getFaciesCategories(float[] lithologies)
    {
        if (lithologies == null) return null;
        float[] facies = new float[lithologies.length];
        for (int i=0; i<lithologies.length; i++)
        {
            Integer lithologyObject = new Integer((int)lithologies[i]);
            String faciesName = (String)lithologyFaciesTable.get(lithologyObject);
            facies[i] = (float)getFaciesCategory(faciesName);
        }
        return facies;
    }

    /** get the most commonly occuring facies for an array of lithologies */
    public float getMostCommonFaciesCategory(float[] lithologies)
    {
        float[] facies = getFaciesCategories(lithologies);
        if (facies == null) return (float)INVALID_FACIES;
        int[] count = new int[nFacies];
        for (int i=0; i<nFacies; i++) count[i] = 0;
        for (int i=0; i<facies.length; i++)
        {
            int intFacies = (int)facies[i];
            if (intFacies == INVALID_FACIES) return INVALID_FACIES;
            count[intFacies]++;
        }
        int maxCount = 0;
        float maxFaciedCategory = (float)INVALID_FACIES;
        for (int i=0; i<nFacies; i++)
        {
            if (count[i] > maxCount)
            {
                maxCount = count[i];
                maxFaciedCategory = i;
            }
        }
        return maxFaciedCategory;
    }

    /** diagnostic output */
    public void print()
    {
        final String NL = System.getProperty("line.separator");

        StringBuffer buffer = new StringBuffer();
        buffer.append("nFacies = " + nFacies + NL);

        // print lithology list
        buffer.append("Lithologies: " + NL);
        int[] lithologies = getAllLithologies();
        if (lithologies == null) buffer.append("< no lithologies >" + NL);
        else
        {
            for (int i=0; i<lithologies.length; i++)
            {
                buffer.append("  " + lithologies[i]);
                if ((i+1)%10 == 0) buffer.append(NL);
            }
        }
        buffer.append(NL);

        // print hash table
        buffer.append("Lithology-Facies Pairs:" + NL);
        int nPairs = lithologyFaciesTable.size();
        if (nPairs < 1) buffer.append("< no lithology-facies pairs >" + NL);
        else
        {
            Enumeration keys = lithologyFaciesTable.keys();
            Enumeration names = lithologyFaciesTable.elements();
            buffer.append("Lithology\tFacies" + NL);
            while (keys.hasMoreElements() && names.hasMoreElements())
            {
                Integer lithologyObject = (Integer)keys.nextElement();
                String faciesName = (String)names.nextElement();
                buffer.append(lithologyObject + "\t" + faciesName + NL);
            }

        }
        buffer.append(NL);

        // print facies lists
        buffer.append("Facies Lithology Lists:" + NL);
        for (int i=0; i<nFacies; i++)
        {
            buffer.append("Facies " + getFaciesName(i) + ":" + NL);
            lithologies = this.getLithologies(i);
            if (lithologies == null)
            {
                buffer.append("< no lithologies >" + NL);
            }
            else
            {
                for (int j=0; j<lithologies.length; j++)
                {
                    buffer.append("  " + lithologies[j]);
                    if ((j+1)%10 == 0) buffer.append(NL);
                }
            }
            buffer.append(NL);
        }
        buffer.append(NL);

        PrintWriter out = new PrintWriter(System.out, true); // needed for correct formatting
        out.println(buffer.toString());
    }


    /** test program */
    public static void main(String[] args)
    {
        final int NFAC = 3;

        final float[] LITH1 = { 5, 11, 3, 57, 22, 5, 5, 3, 5, 1 };
        final float[] LITH2 = { 33, 3, 22, 5, 21, 0, 4, 11 };
        final float[] LITH3 = { 65, 5, 27, 55, 2, 44, 12, 13, 41, 99 };

        final int[] FAC_PR1 =  { 0, 0, 2, 1, 2, 0, 1, 1, 1, 0, 2, 2, 2, 0, 1, 1, 1, 0 };
        final int[] LITH_PR1 = { 5, 3, 21, 0, 3, 57, 11, 22, 1, 33, 4, 65, 27, 55, 2, 44, 12, 13 };

        final int[] LITH_REPLACE = { 99, 2, 5, 1, 41 };

        final int[] FAC_PR2 =  { 0, 0, 1, 2, 2 };
        final int[] LITH_PR2 = { 99, 2, 5, 1, 41 };

        final float[] LITH_VALUES = { 5, 3, 22, 0, 4, 4, 11, 11, 2, 2, 1, 2, 0,
                                    99, 5, 4, 4, 4, 4, 41, 41, 3, 0 };

        final int[] FAC_PR3 = { 1, 2, 0 };
        final int[] LITH_PR3 = { 22, 0, 11 };

        try
        {
            StsCategoricalFacies cf = new StsCategoricalFacies(NFAC);
            cf.addLithologies(LITH1);
            cf.addLithologies(LITH2);
            cf.addLithologies(LITH3);

            for (int i=0; i<FAC_PR1.length; i++)
            {
                cf.set(FAC_PR1[i], LITH_PR1[i]);
            }
            for (int i=0; i<FAC_PR2.length; i++)
            {
                cf.set(FAC_PR2[i], LITH_PR2[i]);
            }
            System.out.println("*** After Pairs Set #2 ***");
            cf.print();

            cf.addReplace(1, LITH_REPLACE);
            System.out.println("\n\n*** After facies B list replaced ***");
            cf.print();

            float[] facies = cf.getFaciesCategories(LITH_VALUES);
            for (int i=0; i<LITH_VALUES.length; i++)
            {
                System.out.println("lith = " + LITH_VALUES[i] +
                        ", facies = " + facies[i] + "\n");
            }

            System.out.println("most common lith = " +
                    cf.getMostCommonFaciesCategory(LITH_VALUES) + "\n");


            for (int i=0; i<FAC_PR3.length; i++)
            {
                cf.set(FAC_PR3[i], LITH_PR3[i]);
            }
            System.out.println("*** After Pairs #3 Set ***");
            cf.print();

            facies = cf.getFaciesCategories(LITH_VALUES);
            for (int i=0; i<LITH_VALUES.length; i++)
            {
                System.out.println("lith = " + (int)LITH_VALUES[i] +
                        ", facies = " + (int)facies[i] + "\n");
            }

            System.out.println("most common lith = " +
                    (int)cf.getMostCommonFaciesCategory(LITH_VALUES) + "\n");

            cf.removeFacies(getFaciesCategory("C"));
            System.out.println("*** After Removing facies C lithologies ***");
            cf.print();
        }
        catch (StsException e) { e.printStackTrace(); }

        System.exit(0);
    }

}




