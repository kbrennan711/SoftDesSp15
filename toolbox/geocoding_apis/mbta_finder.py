"""
Kelly Brennan

Geocoding and Web APIs Project Toolbox exercise

Find the MBTA stops closest to a given location.

Full instructions are at:
https://sites.google.com/site/sd15spring/home/project-toolbox/geocoding-and-web-apis
"""

import urllib   # urlencode function
import urllib2  # urlopen function (better than urllib version)
import json


# Useful URLs (you need to add the appropriate parameters for your requests)
GMAPS_BASE_URL = "https://maps.googleapis.com/maps/api/geocode/json"
MBTA_BASE_URL = "http://realtime.mbta.com/developer/api/v2/stopsbylocation"
MBTA_DEMO_API_KEY = "wX9NwuHnZU2ToO7GmGR9uw"


# A little bit of scaffolding if you want to use it

def get_json(url):
    """
    Given a properly formatted URL for a JSON web API request, return
    a Python JSON object containing the response to that request.
    """
    f = urllib2.urlopen(url)
    response_text = f.read()
    return json.loads(response_text)


def get_lat_long(place_name):
    """
    Given a place name or address, return a (latitude, longitude) tuple
    with the coordinates of the given place.

    See https://developers.google.com/maps/documentation/geocoding/
    for Google Maps Geocode API URL formatting requirements.

    >>> get_lat_long("Fenway Park")
    (42.3466764, -71.0972178)
    """
    text_input = place_name.replace(" ", "")
    response_data = get_json(GMAPS_BASE_URL+"?address=?"+text_input)
    lat = response_data["results"][0]["geometry"]["location"]["lat"]
    lon = response_data["results"][0]["geometry"]["location"]["lng"]
    return (lat, lon)


def get_nearest_station(latitude, longitude):
    """
    Given latitude and longitude strings, return a (station_name, distance)
    tuple for the nearest MBTA station to the given coordinates.

    See http://realtime.mbta.com/Portal/Home/Documents for URL
    formatting requirements for the 'stopsbylocation' API.

    >>> get_nearest_station(str(42.3466764), str(-71.0972178))
    (u'Kenmore', u'0.188732624053955')
    """
    response_data = get_json(MBTA_BASE_URL + "?api_key=" + MBTA_DEMO_API_KEY + "&lat=" + latitude + "&lon=" + longitude)
    for i in range(len(response_data["stop"])): #["distance"]
        # print response_data["stop"][i]["parent_station_name"]
        if response_data["stop"][i]["parent_station_name"] == "":
            pass
        else:
            station_name = response_data["stop"][i]["parent_station_name"]
            distance = response_data["stop"][i]["distance"] 
            return(station_name, distance)
    print "No station within 1 mile"

# print get_nearest_station(str(42.3466764), str(-71.0972178))


def find_stop_near(place_name):
    """
    Given a place name or address, print the nearest MBTA stop and the 
    distance from the given place to that stop.

    >>> find_stop_near("Brigham and Women's Hospital")
    Nearest Station: Brigham Circle
    Distance: 0.128925859928131

    >>> find_stop_near("1000 Olin Way Needham, MA")
    No station within 1 mile
    """
    lat, lon = get_lat_long(place_name)
    if get_nearest_station(str(lat), str(lon)) == None:
        pass
    else:
        station, distance = get_nearest_station(str(lat), str(lon))
        print "Nearest Station:", station
        print "Distance:", distance

find_stop_near("Brigham and Women's Hospital")

if __name__ == '__main__':
    import doctest
    doctest.testmod()
