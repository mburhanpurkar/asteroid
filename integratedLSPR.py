""" Performs centroiding and least squares plate reduction on starfield images.

IMAGE CORRECTIONS: The starfield image is first corrected using a flat image
and the background pixels are subtracted.

CENTROIDING: A weighted mean is taken using the dimensions and estimated centroid
specified by the user.

VISUALIZATION: Quick check centroiding worked properly using vpython - represents
pixel intensity as sphere radius.

ERROR: Uses variance to estimate max and min possible centroid values.

LSPR: Uses least squares plate reduction to find the RA and DEC of the asteroid
given its pixel values and the pixel, RA, and DEC values for six reference stars
in the starfield image.

Lastly, appologies for any poor implementations. This was written back before I
knew how to code well :)

12 July 2015 """


from pyfits import getdata
import matplotlib.pyplot as plt
from numpy import *


""" CENTROIDING """


def showImage(inputImage):
    # Takes as input an object image and displays it
    objectFits = getdata(inputImage)
    plt.imshow(objectFits, interpolation='none')
    plt.gray()
    plt.show()


def adjustedImage(inputImage, inputFlat):
    objectFits = getdata(inputImage)
    # Normalizes the flat image
    normFlat = inputFlat * 1.
    mean1 = normFlat.mean()
    normFlat = normFlat / mean1
    # Divide the object by the normalized flat image
    adjustedObject = objectFits * 1.
    for row in range(len(normFlat)):
        for col in range(len(normFlat[0])):
            if normFlat[row, col] == 0:
                adjustedObject[row, col] = 0
            else:
                adjustedObject[row, col] = objectFits[row, col] / normFlat[row, col]
    return adjustedObject


def darkCompensation(inputImage, adjustedObject):
    objectFits = getdata(inputImage)
    # Slice relevant part of image
    darkArea = objectFits * 1.
    darkArea = darkArea[darkAreay1:darkAreay2, darkAreax1:darkAreax2]
    # Find the mean pixel intensity of the dark area and and subtract it from obj
    meanDark = darkArea.mean()
    for row in range(len(adjustedObject)):
        for col in range(len(adjustedObject[0])):
            adjustedObject[row, col] = adjustedObject[row, col] - meanDark
    return adjustedObject


def findCentroid(adjustedObject, xcoord, ycoord, dimension):
    # User specifies coordinate of pixel, program uses x coordinates to find centroid
    n = (dimension - 1) / 2
    centroidSlice = adjustedObject * 1.
    centroidSlice = centroidSlice[ycoord - n:ycoord + (n + 1), xcoord - n:xcoord + (n + 1)]
    xCenter = len(centroidSlice[0]) / 2
    yCenter = len(centroidSlice) / 2
    weightedSumX = 0
    weightedSumY = 0
    sumOfNumbers = 0
    for row in range(len(centroidSlice)):
        for col in range(len(centroidSlice[0])):
            weightedSumX += centroidSlice[row][col] * (col - xCenter)
            weightedSumY += centroidSlice[row][col] * (yCenter - row)
            sumOfNumbers += centroidSlice[row][col]
    xcentroid = 1. * weightedSumX / sumOfNumbers + xcoord
    ycentroid = -1. * weightedSumY / sumOfNumbers + ycoord
    print "The centroid is " + str(xcentroid) + ", " + str(ycentroid)
    coords = [xcentroid, ycentroid]
    # Old visualization code (to use, from visual import *)
    # x = xcentroid - xcoord
    # y = ycentroid - ycoord
    # image3x3 = copy(centroidSlice)
    # image3x3mean = image3x3.mean()
    # normImage = image3x3 / image3x3mean
    # sphere00 = sphere(pos=(-1,1,0), radius = .1 * normImage[0,0], color = color.blue)
    # sphere01 = sphere(pos=(0,1,0), radius = .1 * normImage[0,1], color = color.blue)
    # sphere02 = sphere(pos=(1,1,0), radius = .1 * normImage[0,2], color = color.blue)
    # sphere10 = sphere(pos=(-1,0,0), radius = .1 * normImage[1,0], color = color.blue)
    # sphere11 = sphere(pos=(0,0,0), radius = .1 * normImage[1,1], color = color.blue)
    # sphere12 = sphere(pos=(1,0,0), radius = .1 * normImage[1,2], color = color.blue)
    # sphere20 = sphere(pos=(-1,-1,0), radius = .1 * normImage[2,0], color = color.blue)
    # sphere21 = sphere(pos=(0,-1,0), radius = .1 * normImage[2,1], color = color.blue)
    # sphere22 = sphere(pos=(1,-1,0), radius = .1 * normImage[2,2], color = color.blue)
    # spherecent = sphere(pos=(x,-1*y,0), radius = .2, color = color.red)
    return coords


""" RUDIMENTARY ERROR ANALYSIS """


def error(coords, img1, img2, img3, img4, img5, xcoord, ycoord, adjustedObject, dimension):
    # Find standard devaition array
    array1 = getdata(img1)
    array2 = getdata(img2)
    array3 = getdata(img3)
    array4 = getdata(img4)
    array5 = getdata(img5)
    listOfArrays = [array1, array2, array3, array4, array5]
    averageArray = zeros(shape(array1))
    varianceArray = zeros(shape(array1))
    for a in listOfArrays:
        averageArray += a
    averageArray = averageArray / 5.
    for a in listOfArrays:
        varianceArray += square(a - averageArray)
    varianceArray = sqrt(varianceArray / 5.)
    # Find x max and min data array
    n = (dimension - 1) / 2
    smallVariance = varianceArray * 1.
    smallVariance = smallVariance[ycoord - n:ycoord + (n + 1), xcoord - n:xcoord + (n + 1)]
    smallDatax = adjustedObject * 1.
    smallDatax = adjustedObject[ycoord - n:ycoord + (n + 1), xcoord - n:xcoord + (n + 1)]
    centerx = len(smallVariance[0]) / 2
    centery = len(smallVariance) / 2
    for row in range(len(smallVariance)):
        for col in range(len(smallVariance[0])):
            if col < centerx:
                smallDatax[row][col] -= smallVariance[row][col]
            elif col > centerx:
                smallDatax[row][col] += smallVariance[row][col]
    # smallDatax is now the x max and min array
    # Now, find the centroid of this new array
    weightedSumX = 0
    weightedSumY = 0
    sumOfNumbers = 0
    for row in range(len(smallDatax)):
        for col in range(len(smallDatax[0])):
            weightedSumX += smallDatax[row][col] * (col - centerx)
            sumOfNumbers += smallDatax[row][col]
    xcentroid = abs((1. * weightedSumX / sumOfNumbers + xcoord) - coords[0])
    print "The x uncertainty is " + str(xcentroid)
    # Find y max and min data array
    smallDatay = adjustedObject * 1.
    smallDatay = adjustedObject[ycoord - n:ycoord + (n + 1), xcoord - n:xcoord + (n + 1)]
    for row in range(len(smallVariance)):
        for col in range(len(smallVariance[0])):
            if row < centery:
                smallDatay[row][col] += smallVariance[row][col]
            elif row > centery:
                smallDatay[row][col] -= smallVariance[row][col]
    # smallDatay is now the y max and min array
    # Now, find the centroid of this new array
    weightedSumX = 0
    weightedSumY = 0
    sumOfNumbers = 0
    for row in range(len(smallDatay)):
        for col in range(len(smallDatay[0])):
            weightedSumY += smallDatay[row][col] * (centery - row)
            sumOfNumbers += smallDatay[row][col]
    ycentroid = abs((1. * weightedSumX / sumOfNumbers + ycoord) - coords[1])
    print "The y uncertainty is " + str(ycentroid)
    return xcentroid, ycentroid


""" LSPR """


def intensity(first_image, starx, stary, starBuffer, blankx, blanky):
    # This calculates the pixel counts of stars
    star = getdata(first_image)
    starbackground = star * 1.
    starbackground = starbackground[(stary - starBuffer):(stary + starBuffer), (starx - starBuffer):(starx + starBuffer)]
    starbackground = sum(starbackground)
    background = star * 1.
    background = background[(blanky - starBuffer):(blanky + starBuffer), (blankx - starBuffer):(blankx + starBuffer)]
    background = sum(background)
    brightObject = starbackground - background
    return brightObject


def determinant2x2(arr):
    # Inputs a two by two array and outputs the determinant
    a = arr[0][0]
    b = arr[0][1]
    c = arr[1][0]
    d = arr[1][1]
    determinant = a * d - b * c
    return float(determinant)


def determinant3x3(arr):
    # Inputs a three by three array and outputs the determinant
    determinant = 0
    for i in range(0, 3):
        newarr = append(arr[1:, 0:i], arr[1:, i + 1:], axis=1)
        determinant += ((-1)**i) * arr[0][i] * determinant2x2(newarr)
    return float(determinant)


def cramer3x3(coeff, const):
    # Inputs a three by three array and a one by three array and outputs the solutions
    coeffdeter = float(determinant3x3(coeff))
    xarr = array(coeff)
    xarr[:, 0] = const
    xdeter = determinant3x3(xarr)
    yarr = array(coeff)
    yarr[:, 1] = const
    ydeter = determinant3x3(yarr)
    zarr = array(coeff)
    zarr[:, 2] = const
    zdeter = determinant3x3(zarr)
    x = float(xdeter) / coeffdeter
    y = float(ydeter) / coeffdeter
    z = float(zdeter) / coeffdeter
    deterList = [x, y, z]
    return deterList


def raDecCalculator(arr, asteroid, firstImage, n, starXY):
    # This function calculates the sum of quantities in a coeffecient array to
    # prepare for a regression calculation. It then does the regression.
    # Note that arr takes in values that are
    # in the form (RA, DEC, x, y). It then creates two coeffecient arrays
    x, x2, xy, y, y2, ax, ay, a, bx, by, b = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # Does the summation calculations
    for i in range(len(arr)):
        x += arr[i, 2]
        x2 += (arr[i, 2])**2
        xy += (arr[i, 2]) * (arr[i, 3])
        y += arr[i, 3]
        y2 += (arr[i, 3])**2
        ax += (arr[i, 0]) * (arr[i, 2])
        ay += (arr[i, 0]) * (arr[i, 3])
        a += arr[i, 0]
        bx += (arr[i, 1]) * (arr[i, 2])
        by += (arr[i, 1]) * (arr[i, 3])
        b += arr[i, 1]
    # Makes the coeffecient array
    coefficientArray = array([[len(arr), x, y],
                              [x, x2, xy],
                              [y, xy, y2]])
    # Constant RA array
    smallerRAArray = array([[a, ax, ay]])
    # RA solutions
    solutionsRA = cramer3x3(coefficientArray, smallerRAArray)
    # Constant DEC array
    smallerDECArray = array([[b, bx, by]])
    # DEC Solutions
    solutionsDEC = cramer3x3(coefficientArray, smallerDECArray)
    # Calculates values of b1, a11, a12, b2, a21, a22
    b1 = solutionsRA[0]
    a11 = solutionsRA[1]
    a12 = solutionsRA[2]
    b2 = solutionsDEC[0]
    a21 = solutionsDEC[1]
    a22 = solutionsDEC[2]
    # Find RA and DEC
    ra = b1 + a11 * asteroid[0] + a12 * asteroid[1]
    dec = b2 + a21 * asteroid[0] + a22 * asteroid[1]
    print "\nAsteroid RA: " + str(ra)
    print "Asteroid Dec: " + str(dec)
    # Find residuals
    acceptedRA = (arr[:, 0]).tolist()
    acceptedDEC = (arr[:, 1]).tolist()
    calculatedRA = []
    calculatedDEC = []
    residsRA = []
    residsDEC = []
    for i in range(len(arr)):
        ra = 0
        dec = 0
        ra += (b1 + a11 * arr[i, 2] + a12 * arr[i, 3])
        dec += (b2 + a21 * arr[i, 2] + a22 * arr[i, 3])
        calculatedRA.append(ra)
        calculatedDEC.append(dec)
    residsRA = []
    residsDEC = []
    for i in range(len(arr)):
        rra = 0
        rdec = 0
        rra = (acceptedRA[i] - calculatedRA[i])
        rdec = (acceptedDEC[i] - calculatedDEC[i])
        residsRA.append(rra)
        residsDEC.append(rdec)
    # Do variance stuff
    residsRAsum = 0
    residsDECsum = 0
    for i in range(len(arr)):
        residsRAsum += (residsRA[i])**2
        residsDECsum += (residsDEC[i])**2
    raSigma = sqrt((residsRAsum) / 3)
    decSigma = sqrt((residsDECsum) / 3)
    raVariance = residsRAsum / 3
    decVariance = residsDECsum / 3
    print "\nRA Residuals"
    print residsRA
    print "\nDEC Residuals"
    print residsDEC
    print "\nRA Variance: " + str(raVariance)
    print "DEC Variance: " + str(decVariance)
    print "RA Standard Deviation: " + str(raSigma)
    print "DEC Standard Deviation: " + str(decSigma)
    # Find the intensity of stars and asteroid - loop through all 6-8 stars
    # Creates a list with 6-8 star intensities to be used as proportions in visualization
    blankx = input("Enter the x coordinate of a blank area: ")
    blanky = input("Enter the y coordinate of a blank area: ")
    intensityList = []
    for element in range(len(starXY)):
        starx = starXY[element, 2]
        stary = starXY[element, 3]
        starBuffer = input("Enter a buffer around the " + str(element + 1) + " star: ")
        sIntensity = intensity(firstImage, starx, stary, starBuffer, blankx, blanky)
        intensityList += sIntensity
    starsValuesRadec = []
    starsValuesRa = []
    starsValuesDec = []
    for i in range(len(arr)):
        starsValuesRadec.append(arr[i, 0])
        starsValuesRadec.append(arr[i, 1])
    for i in range(len(arr)):
        starsValuesRa.append(starsValuesRadec[2 * i])
        starsValuesDec.append(starsValuesRadec[2 * i + 1])
    for i in range(len(intensityList)):
        sphere(pos=(starsValuesRa[i], starsValuesDec[i], 0), radius=.05 * intensityList[i])
    asteroid = box(pos=(ra, dec, 0), length=.05 * raSigma, height=.05 * decSigma)
    return ra, dec


""" INPUTS """

# User inputs images
flatInput = raw_input("Enter the flat image's file name: ")
flat = getdata(flatInput)
firstImage = raw_input("Enter the first image's file name: ")
secondImage = raw_input("Enter the second image's file name: ")
thirdImage = raw_input("Enter the third image's file name: ")
fourthImage = raw_input("Enter the fourth image's file name: ")
fifthImage = raw_input("Enter the fifth image's file name: ")
# Image is shown
showImage(firstImage)
# User inputs things for calculating centroids of stars
i = input("How many stars would you like to use in your calculation? ")
n = i
darkAreax1 = input("Enter the first x coordinate of the dark area: ")
darkAreax2 = input("Enter the second x coordinate of the dark area: ")
darkAreay1 = input("Enter the first y coordinate of the dark area: ")
darkAreay2 = input("Enter the second y coordinate of the dark area: ")
dimension = input("What is the length of your centroid array? ")
# User enters 6 to 8 stars, program calculates centroid 6 to 8 times
# These values are imported into the LSPR program, which takes in the coordinates
# of the asteroid in x and y, as well as the RA and DEC of the stars, and uses
# them to calculate the asteroid centroid
outputArray = []
while i > 0:
    brightX = input('What is the x coordinate of a star? ')
    brightY = input('What is the y coordinate of a star? ')
    RA = input("What is the RA of the star? ")
    dec = input("What is the DEC of the star? ")
    adjusted = adjustedImage(firstImage, flat)
    compensated = darkCompensation(firstImage, adjusted)
    centroid = findCentroid(compensated, brightX, brightY, dimension)
    RA = [RA]
    dec = [dec]
    temp = RA + dec + centroid
    outputArray += temp
    error(centroid, firstImage, secondImage, thirdImage, fourthImage, fifthImage, brightX, brightY, compensated, dimension)
    i -= 1
# Array with star RA, dec, x, and y values
outputArray = array(outputArray)
outputArray.shape = ((len(outputArray) / 4), 4)
# User inputs asteroid x and y coordinates, and the centroiding program is run on them
# to find the real centroid of the asteroid. The asteroid RA and DEC are then calculated
asteroid = input("Please input the coordinates of your asteroid (as a list): ")
centroidAsteroid = findCentroid(compensated, asteroid[0], asteroid[1], dimension)
asteroidCoords = raDecCalculator(outputArray, centroidAsteroid, firstImage, n, outputArray)
