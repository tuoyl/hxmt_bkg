#!/usr/bin/env python
 

print "*********** Detector ID catalogue ***********"

payloads = raw_input("what do you want?(HE/ME/LE):\n")

if payloads == "HE" or payloads == "he":
    print "\nHE small detector ID:"
    print "0, 1, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17"
    print "--------------------------"
    print "\nHE blind detector ID:"
    print "16"
    print "--------------------------"

   
elif payloads == "ME" or payloads == "me":
    print "\nME small detector ID:"
    print "0-5, 7, 12-23, 25, 30-41, 43, 48-53"
    print "--------------------------"
    print "\nME blind detector ID:"
    print "10, 28, 46"
    print "--------------------------"
elif payloads == "LE" or payloads == "le":
    print "\nLE small detector ID:"
    print '0,2-4,6-10,12,14,20,22-26,28-30,32,34-36,'+\
            '38-42,44,46,52,54-58,60-62'+\
            '64,66-68,70-74,76,78,84,86-90,92-94'
    print "--------------------------"
    print "\nLE blind detector ID:"
    print "13, 45, 77"
    print "--------------------------"

