package com.Sts.Utilities.DateTime;

import junit.framework.*;

import java.util.*;

class TinyDate
{
    private static final String[] MONTHNAME =
            {
                    "January", "February", "March", "April", "May", "June",
                    "July", "August", "September", "October", "November", "December"
            };

    private String str;
    private int order, year, month, date, hour, minute, second, milli;

    TinyDate(String str, int order, int year, int month, int date)
    {
        this(str, order, year, month, date, 0, 0, 0, 0);
    }

    TinyDate(String str, int order, int year, int month, int date, int hour,
             int minute, int second)
    {
        this(str, order, year, month, date, hour, minute, second, 0);
    }

    TinyDate(String str, int order, int year, int month, int date, int hour,
             int minute, int second, int milli)
    {
        this.str = str;
        this.order = order;
        this.year = year;
        this.month = month - 1;
        this.date = date;
        this.hour = hour;
        this.minute = minute;
        this.second = second;
        this.milli = milli;
    }

    public boolean checkDate(Calendar cal)
    {
        return (year == cal.get(Calendar.YEAR) &&
                month == cal.get(Calendar.MONTH) &&
                date == cal.get(Calendar.DATE) &&
                hour == cal.get(Calendar.HOUR_OF_DAY) &&
                minute == cal.get(Calendar.MINUTE) &&
                second == cal.get(Calendar.SECOND) &&
                milli == cal.get(Calendar.MILLISECOND));
/*
if (year == cal.get(Calendar.YEAR) && month == cal.get(Calendar.MONTH) &&
    date == cal.get(Calendar.DATE) && hour == cal.get(Calendar.HOUR_OF_DAY) &&
    minute == cal.get(Calendar.MINUTE) && second == cal.get(Calendar.SECOND) &&
    milli == cal.get(Calendar.MILLISECOND))
{
    return true;
}
System.err.println("CHKDATE: " + cal.get(Calendar.YEAR) + "/" + cal.get(Calendar.MONTH) + "/" + cal.get(Calendar.DATE) + " " + cal.get(Calendar.HOUR_OF_DAY) + ":" + cal.get(Calendar.MINUTE) + ":" + cal.get(Calendar.SECOND) + "." + cal.get(Calendar.MILLISECOND) + " (" + CalendarParser.toString(cal) + ") != " + year + "/" + month + "/" + date + " " + hour + ":" + minute + ":" + second + "." + milli + " (" + toString() + ")");
return false;
*/
    }

    public int getOrder()
    {
        return order;
    }

    public String getString()
    {
        return str;
    }

    public String toString()
    {
        final String timeStr;
        if(milli > 0)
        {
            timeStr = " " + hour +
                    (minute < 10 ? ":0" : ":") + minute +
                    (second < 10 ? ":0" : ":") + second +
                    (milli < 100 ? ".00" : (milli < 10 ? ".0" : ".")) + milli;
        }
        else if(second > 0)
        {
            timeStr = " " + hour +
                    (minute < 10 ? ":0" : ":") + minute +
                    (second < 10 ? ":0" : ":") + second;
        }
        else if(hour > 0 || minute > 0)
        {
            timeStr = " " + hour +
                    (minute < 10 ? ":0" : ":") + minute;
        }
        else
        {
            timeStr = "";
        }

        final String zoneStr;
        if(timeStr.length() == 0)
        {
            zoneStr = "";
        }
        else
        {
            zoneStr = " GMT";
        }

        return "" + MONTHNAME[month] + " " + date + ", " + year + timeStr +
                zoneStr;
    }
}

public class CalendarParserTest extends TestCase
{
    private boolean fullTest = (System.getProperty("fullTest") != null);

    public CalendarParserTest(String name)
    {
        super(name);
    }

    private static Calendar checkAll(String testStr, int order, boolean isValid)
    {
        return checkAll(testStr, order, isValid, false);
    }

    private static Calendar checkAll(String testStr, int order, boolean isValid, boolean ignoreChanges)
    {
        Calendar cal;
        try
        {
            cal = CalendarParser.parse(testStr, order, ignoreChanges);
        }
        catch(Throwable e)
        {
            if(isValid)
            {
                fail("\"" + testStr + "\" (order " +
                        CalendarParser.getOrderString(order) + ") threw " +
                        e.getClass().getName() + ": " + e.getMessage());
            }
            return null;
        }

        if(!isValid)
        {
            fail("Didn't expect \"" + testStr + "\" to make date " + cal);
        }

        String str = CalendarParser.toString(cal);

        Calendar cal2;
        try
        {
            cal2 = CalendarParser.parse(str);
        }
        catch(Throwable e)
        {
            if(isValid)
            {
                fail("\"" + str + "\" (toString \"" + testStr + "\") threw " +
                        e.getClass().getName() + ": " + e.getMessage());
            }
            cal2 = null;
        }

        assertTrue("String \"" + str + "\" doesn't match \"" +
                CalendarParser.toString(cal2), equals(cal, cal2));

        String pretty = CalendarParser.prettyString(cal);
        try
        {
            cal2 = CalendarParser.parse(pretty);
        }
        catch(Throwable e)
        {
            if(isValid)
            {
                fail("\"" + pretty + "\" (prettyString \"" + testStr +
                        "\") threw " + e.getClass().getName() + ": " +
                        e.getMessage());
            }
            cal2 = null;
        }

        assertTrue("Pretty \"" + pretty + "\" doesn't match \"" +
                CalendarParser.prettyString(cal2), equals(cal, cal2));

        String sql = CalendarParser.toSQLString(cal);
        try
        {
            cal2 = CalendarParser.parse(sql);
        }
        catch(Throwable e)
        {
            if(isValid)
            {
                fail("\"" + sql + "\" (sqlString \"" + testStr + "\") threw " +
                        e.getClass().getName() + ": " + e.getMessage());
            }
            cal2 = null;
        }

        assertTrue("SQL \"" + sql + "\" doesn't match \"" +
                CalendarParser.toSQLString(cal2), equals(cal, cal2));

        return cal2;
    }

    private static final int compare(Calendar c1, Calendar c2)
    {
        if(c1 == null)
        {
            return (c2 == null ? 0 : 1);
        }
        else if(c2 == null)
        {
            return -1;
        }

        int tCmp = c2.get(Calendar.YEAR) - c1.get(Calendar.YEAR);
        if(tCmp != 0)
        {
            tCmp = c2.get(Calendar.MONTH) - c1.get(Calendar.MONTH);
            if(tCmp != 0)
            {
                tCmp = c2.get(Calendar.DATE) - c1.get(Calendar.DATE);
                if(tCmp != 0)
                {
                    tCmp = c2.get(Calendar.HOUR) - c1.get(Calendar.HOUR);
                    if(tCmp != 0)
                    {
                        tCmp = c2.get(Calendar.MINUTE) -
                                c1.get(Calendar.MINUTE);
                        if(tCmp != 0)
                        {
                            tCmp = c2.get(Calendar.SECOND) -
                                    c1.get(Calendar.SECOND);
                            if(tCmp != 0)
                            {
                                tCmp = c2.get(Calendar.MILLISECOND) -
                                        c1.get(Calendar.MILLISECOND);
                            }
                        }
                    }
                }
            }
        }

        if(tCmp != 0)
        {
            return (tCmp < 0 ? -1 : 1);
        }

        TimeZone z1 = c1.getTimeZone();
        TimeZone z2 = c2.getTimeZone();

        final int zCmp = z2.getRawOffset() - z1.getRawOffset();
        if(zCmp != 0)
        {
            return (zCmp < 0 ? -1 : 1);
        }

        return 0;
    }

    private static final boolean equals(Calendar c1, Calendar c2)
    {
        return (compare(c1, c2) == 0);
    }

    private static final String getTimeZoneString()
    {
        TimeZone tz = TimeZone.getDefault();

        int raw = tz.getRawOffset() / (60 * 1000);

        StringBuffer buf = new StringBuffer();
        if(raw >= 0)
        {
            buf.append('+');
        }
        else
        {
            buf.append('-');
            raw = -raw;
        }

        final int hour = raw / 60;
        if(hour < 10)
        {
            buf.append('0');
        }
        buf.append(hour);

        buf.append(':');

        final int min = raw % 60;
        if(min < 10)
        {
            buf.append('0');
        }
        buf.append(min);

        return buf.toString();
    }

    public static Test suite()
    {
        return new TestSuite(CalendarParserTest.class);
    }

    public void testNull()
            throws CalendarParserException
    {
        assertNull("Didn't expect CalendarParser.parse(null) to return anything!",
                checkAll(null, CalendarParser.DD_MM_YY, true));
    }

    public void testPermutations()
    {
        if(!fullTest)
        {
            return;
        }

        final String[] monthName =
                {
                        "January",
                        "February",
                        "March",
                        "April",
                        "May",
                        "June",
                        "July",
                        "August",
                        "September",
                        "October",
                        "November",
                        "December"
                };

        for(int m = 0; m < 37; m++)
        {

            final int mNum;
            final String mStr;
            final boolean mValid;
            if(m < 12)
            {
                mNum = m + 1;
                mStr = monthName[m];
                mValid = true;
            }
            else if(m < 24)
            {
                mNum = (m - 12) + 1;
                mStr = monthName[mNum - 1].substring(0, 3);
                mValid = true;
            }
            else if(m < 36)
            {
                mNum = (m - 24) + 1;
                mStr = Integer.toString(mNum);
                mValid = true;
            }
            else
            {
                mNum = -1;
                mStr = "Barf";
                mValid = false;
            }

            for(int d = 1; d < 33; d++)
            {
                final int dNum = d;
                final String dStr = Integer.toString(dNum);
                final boolean dValid;
                if(d < 1 || d > 31)
                {
                    dValid = false;
                }
                else if(d > 28 && mNum == 2)
                {
                    dValid = false;
                }
                else if(d > 30 && (mNum == 4 || mNum == 6 || mNum == 9 ||
                        mNum == 11))
                {
                    dValid = false;
                }
                else
                {
                    dValid = true;
                }

                for(int y = 0; y < 4; y++)
                {
                    final int yNum;
                    final String yStr;
                    final boolean yValid;
                    switch(y)
                    {
                        case 0:
                            yNum = 2001;
                            yStr = "2001";
                            yValid = true;
                            break;
                        case 1:
                            yNum = 1999;
                            yStr = "1999";
                            yValid = true;
                            break;
                        case 2:
                            yNum = 2001;
                            yStr = "01";
                            yValid = true;
                            break;
                        default:
                            yNum = 1999;
                            yStr = "99";
                            yValid = true;
                            break;
                    }

                    for(int s = 0; s < 3; s++)
                    {
                        final String sepChar;
                        final boolean sValid;
                        switch(s)
                        {
                            case 0:
                                sepChar = " ";
                                sValid = true;
                                break;
                            case 1:
                                sepChar = "-";
                                sValid = true;
                                break;
                            case 2:
                                sepChar = "/";
                                sValid = true;
                                break;
                            default:
                                sepChar = "*";
                                sValid = false;
                                break;
                        }

                        for(int o = 0; o < 8; o++)
                        {
                            final String dateStr;
                            switch(o)
                            {
                                case CalendarParser.DD_MM_YY:
                                    dateStr = dStr + sepChar + mStr + sepChar +
                                            yStr;
                                    break;
                                case CalendarParser.DD_YY_MM:
                                    dateStr = dStr + sepChar + yStr + sepChar +
                                            mStr;
                                    break;
                                case CalendarParser.MM_DD_YY:
                                    dateStr = mStr + sepChar + dStr + sepChar +
                                            yStr;
                                    break;
                                case CalendarParser.MM_YY_DD:
                                    dateStr = mStr + sepChar + yStr + sepChar +
                                            dStr;
                                    break;
                                case CalendarParser.YY_DD_MM:
                                    dateStr = yStr + sepChar + dStr + sepChar +
                                            mStr;
                                    break;
                                case CalendarParser.YY_MM_DD:
                                    dateStr = yStr + sepChar + mStr + sepChar +
                                            dStr;
                                    break;
                                default:
                                    dateStr = null;
                                    break;
                            }

                            if(dateStr != null)
                            {
                                TinyDate td =
                                        new TinyDate(dateStr, o, yNum, mNum, dNum);

                                final boolean dKludge =
                                        (!dValid && yStr.equals("01") && d > 31);

                                final boolean isValid =
                                        mValid && (dValid || dKludge) &&
                                                yValid && sValid;

                                Calendar cal = checkAll(dateStr, o, isValid);
                                if(cal == null)
                                {
                                    assertTrue("Expected \"" + dateStr +
                                            "\" (order " +
                                            CalendarParser.getOrderString(o) +
                                            ") to succeed!",
                                            !isValid);
                                }
                                else if(!isValid)
                                {
                                    fail("Didn't expect \"" + dateStr +
                                            "\" (order " +
                                            CalendarParser.getOrderString(o) +
                                            ") to return \"" +
                                            CalendarParser.toString(cal) +
                                            "\" (bad" +
                                            (mValid ? "" : " month") +
                                            (dValid ? "" : " day") +
                                            (yValid ? "" : " year") +
                                            (sValid ? "" : " separator") +
                                            ")");
                                }
                                else if(!dKludge)
                                {
                                    assertTrue("String \"" + dateStr +
                                            "\" (order " +
                                            CalendarParser.getOrderString(o) +
                                            ") returned " +
                                            cal.get(Calendar.YEAR) + "/" +
                                            cal.get(Calendar.MONTH) + "/" +
                                            cal.get(Calendar.DATE) + " (" +
                                            CalendarParser.toString(cal) +
                                            "), not " + td,
                                            td.checkDate(cal));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void testStrings()
            throws CalendarParserException
    {
        final String tzStr = " " + getTimeZoneString();

        for(int i = 1; i <= 12; i *= 10)
        {
            final String rawStr = "" + i + "-" + i + "-" + i + " " +
                    i + ":" + i + ":" + i + (i == 1 ? ".001" : ".100");


            Calendar cal = CalendarParser.parse(rawStr,
                    CalendarParser.YY_MM_DD,
                    false);

            final String expStd, expPretty, expSQL;
            if(i == 1)
            {
                expStd = "1 Jan 2001  1:01:01.001" + tzStr;
                expPretty = "January 1, 2001  1:01:01.001" + tzStr;
                expSQL = "2001-01-01  1:01:01.001" + tzStr;
            }
            else
            {
                expStd = "10 Oct 2010 10:10:10.100" + tzStr;
                expPretty = "October 10, 2010 10:10:10.100" + tzStr;
                expSQL = "2010-10-10 10:10:10.100" + tzStr;
            }

            assertEquals(rawStr + " standard string mismatch",
                    expStd, CalendarParser.toString(cal));
            assertEquals(rawStr + " pretty string mismatch",
                    expPretty, CalendarParser.prettyString(cal));
            assertEquals(rawStr + " SQL string mismatch",
                    expSQL, CalendarParser.toSQLString(cal));
        }
    }

    public void testBad()
    {
        final String[] badArgs = {
                "",
                "                           ",
                "(Some comment)",
                "Feb",
                "22",
                "1999",
                // ISO 8601 formats
                "2000-0830",
                "2000-830",
                "200008-30",
                "200008",
                "2000-8",
                "2000",
                "00-0830",
                "00-08",
                "00",
                "2000W0202",
                "2000-W02-2",
                "00-W0202",
                "2000045",
                "00-45",
                "2000-08-30T12:34:56Z",
                "2000-08-30T12:34:56+01:00",
                "2000-08-30T12:34:56-01:00",
                // misc formats
                "today",
                "1st thursday in June 1992",
                "december tenth",
                "sunday week 22 1995",
                "22nd sunday",
                "sunday 22nd week in 1996",
                "last friday",
                "next friday",
                "last month",
                "next month",
                "last month in 2000",
                "next month in 2000",
                "3 weeks ago",
                "Friday 2 weeks ago",
                "2 weeks ago friday",
                "last day in October",
                "first day of October",
                "friday",
                // bogus times
                "Mar 3, 2000 1:32 BM",
                "May 8, 2000 1:28 PM RXT",
        };

        for(int i = 0; i < badArgs.length; i++)
        {
            assertNull("Didn't expect \"" + badArgs[i] + "\" to succeed!",
                    checkAll(badArgs[i], CalendarParser.YY_MM_DD, false));
        }

        final String[] moreBad = new String[]{
                "Feb 22",
                "2005 22",
                "2005 Feb",
        };

        for(int i = 0; i < moreBad.length; i++)
        {
            assertNull("Didn't expect \"" + moreBad[i] + "\" to succeed!",
                    checkAll(moreBad[i], CalendarParser.DD_MM_YY, false));
        }
    }

    public void testGood()
    {
        final TinyDate[] goodArgs = {
                new TinyDate("Feb 3 1000", CalendarParser.MM_DD_YY, 1000, 2, 3),
                new TinyDate("February 3, 2000",
                        CalendarParser.MM_DD_YY, 2000, 2, 3),
                new TinyDate("3 Feb 2000", CalendarParser.MM_DD_YY, 2000, 2, 3),
                new TinyDate("2000 Feb 3", CalendarParser.MM_DD_YY, 2000, 2, 3),
                new TinyDate("July 4, 1999", CalendarParser.MM_DD_YY, 1999, 7, 4),
                new TinyDate("Tue, Jan 05, 1999",
                        CalendarParser.MM_DD_YY, 1999, 1, 5),
                new TinyDate("2000/05/08", CalendarParser.MM_DD_YY, 2000, 5, 8),
                new TinyDate("2-3-01", CalendarParser.MM_DD_YY, 2001, 2, 3),
                new TinyDate("1-1-2001", CalendarParser.MM_DD_YY, 2001, 1, 1),
                new TinyDate("2001-2-1", CalendarParser.YY_DD_MM, 2001, 1, 2),
                new TinyDate("Feb 1st, 2001", CalendarParser.MM_DD_YY, 2001, 2, 1),
                new TinyDate("Feb 3rd, 2001", CalendarParser.MM_DD_YY, 2001, 2, 3),
                new TinyDate("Dec 5th, 2001",
                        CalendarParser.MM_DD_YY, 2001, 12, 5),
                new TinyDate("Oct 12th, 2001",
                        CalendarParser.MM_DD_YY, 2001, 10, 12),
                new TinyDate("Oct 12, 2001",
                        CalendarParser.MM_DD_YY, 2001, 10, 12),
                new TinyDate("2001 Oct 12", CalendarParser.MM_DD_YY, 2001, 10, 12),
                new TinyDate("12 2001 Oct", CalendarParser.MM_DD_YY, 2001, 10, 12),
                new TinyDate("1 12 Oct", CalendarParser.MM_DD_YY, 2001, 10, 12),
                new TinyDate("1 12 Oct", CalendarParser.YY_MM_DD, 2001, 10, 12),
                new TinyDate("12/27/00", CalendarParser.MM_DD_YY, 2000, 12, 27),
                new TinyDate("8/30/00", CalendarParser.MM_DD_YY, 2000, 8, 30),
                new TinyDate("8/30/99", CalendarParser.MM_DD_YY, 1999, 8, 30),
                new TinyDate("8/30/49", CalendarParser.MM_DD_YY, 2049, 8, 30),
                new TinyDate("8/30/50", CalendarParser.MM_DD_YY, 1950, 8, 30),
                // ISO 8601 formats
                new TinyDate("2000-08-30", CalendarParser.YY_MM_DD, 2000, 8, 30),
                new TinyDate("2000-8-30", CalendarParser.YY_MM_DD, 2000, 8, 30),
                new TinyDate("20000830", CalendarParser.YY_MM_DD, 2000, 8, 30),
                // times
                new TinyDate("2000-05-08 13:28:19.654",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 13, 28, 19, 654),
                new TinyDate("2000-05-08 1:28:19 PM",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 13, 28, 19),
                new TinyDate("2000-05-08 1:28PM",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 13, 28, 0),
                new TinyDate("2000-05-08 1:28:19 AM",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 1, 28, 19),
                new TinyDate("2000-05-08 1:28:19 AM",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 1, 28, 19),
                new TinyDate("May 8, 2000 1:28 PM CST",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 13, 28, 0),
                new TinyDate("May 8, 2000 7:28 PM GMT",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 19, 28, 0),
                new TinyDate("May 8, 2000 13:28:57 -06:00",
                        CalendarParser.MM_DD_YY, 2000, 5, 8, 13, 28, 57),
        };

        for(int i = 0; i < goodArgs.length; i++)
        {
            Calendar cal = checkAll(goodArgs[i].getString(),
                    goodArgs[i].getOrder(), true);
            assertNotNull("Didn't expect \"" + goodArgs[i] + "\" to fail!",
                    cal);
            assertTrue("String \"" + goodArgs[i].getString() + "\" returned " +
                    cal.get(Calendar.YEAR) + "/" +
                    cal.get(Calendar.MONTH) + "/" +
                    cal.get(Calendar.DATE) + " (" +
                    CalendarParser.toString(cal) + "), not " +
                    goodArgs[i], goodArgs[i].checkDate(cal));
        }
    }

    public void testChanges()
    {
        TinyDate td = new TinyDate("2/30/00",
                CalendarParser.MM_DD_YY, 2000, 3, 1);
        Calendar cal = checkAll(td.getString(), td.getOrder(), true, true);
        assertNotNull("Didn't expect \"" + td + "\" to fail!", cal);
        assertTrue("String \"" + td.getString() + "\" returned " +
                cal.get(Calendar.YEAR) + "/" +
                cal.get(Calendar.MONTH) + "/" +
                cal.get(Calendar.DATE) + " (" +
                CalendarParser.toString(cal) + "), not " +
                td, td.checkDate(cal));
    }

    public static final void main(String[] args)
    {
        if(args.length == 0)
        {
            junit.textui.TestRunner.run(suite());
        }
        else
        {

            boolean result = true;

            for(int i = 0; i < args.length; i++)
            {
                Calendar cal = checkAll(args[i], CalendarParser.YY_MM_DD, true);
                System.err.println(args[i] + " => " + CalendarParser.prettyString(cal));
                long time = cal.getTimeInMillis();
                System.err.println("time: " + time);
                result &= (cal != null);
            }

            if(!result)
            {
                System.exit(1);
            }

            System.out.println("All tests succeeded!");
        }
    }
}
