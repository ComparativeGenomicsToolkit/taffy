// test function for increment, bulk operations

//Database functions
#include <ktremotedb.h>
#include <kclangc.h>

using namespace std;
using namespace kyototycoon;

// main routine
int main(int argc, char** argv) {

  // create the database object
  RemoteDB *rdb = new RemoteDB();

  // open the database
  if (!rdb->open("kolossus-10.kilokluster.ucsc.edu", 1978, -1)) {
    cerr << "open error: " << rdb->error().name() << endl;
  }

  rdb->clear();
  string value;
  // arbitrary expiration time
  int64_t xt = 50;

  // test integer increment function
  int64_t key = 7;
  size_t sizeOfKey = sizeof(uint64_t);
  int64_t initialValue = 123;
  int64_t incrValue = 23;
  int64_t sizeOfRecord = sizeof(uint64_t);
  size_t sp;

  char *record;
  int64_t currValue = 0;
 
  cout << "adding record " << initialValue << " of size " << sizeOfRecord << endl;


  // Normalize a 64-bit number in the native order into the network byte order.
  // little endian (our x86 linux machine) to big Endian....
  int64_t KCSafeIV = kyotocabinet::hton64(initialValue);

  rdb->add((char *)&key, sizeOfKey, (const char *)&KCSafeIV, sizeOfRecord);
  cerr << "add record error: " << rdb->error().name() << endl;

  record = rdb->get((char *)&key, sizeOfKey, &sp, NULL);
  currValue = kyotocabinet::ntoh64(*((int64_t*)record));

  cout << "added record " << currValue << " of size " << sp << endl;

  rdb->increment((char *)&key, sizeOfKey, incrValue, xt);
  rdb->increment((char *)&key, sizeOfKey, incrValue, xt);
  //450 (the existing record was not compatible).
  cerr << "increment error: " << rdb->error().name() << endl;
  record = rdb->get((char *)&key, sizeOfKey, &sp, NULL);

  // Denormalize a 64-bit number in the network byte order into the native order.
  // (big-endian to little endian)
  currValue = kyotocabinet::ntoh64(*((int64_t*)record));

  printf("value after increment (should be 146): %d and size %d\n", currValue, sp);

  rdb->clear();


  printf("trying a second version of increment (flip the byte order of the incr value instead..\n");
  // increment 2 hack:
  rdb->add((char *)&key, sizeOfKey, (const char *)&initialValue, sizeOfRecord);
  cerr << "add record error: " << rdb->error().name() << endl;

  record = rdb->get((char *)&key, sizeOfKey, &sp, NULL);
  currValue = *((uint64_t*)record);

  cout << "added record " << currValue << " of size " << sp << endl;

  int retVal = rdb->increment((char *)&key, sizeOfKey, kyotocabinet::hton64(incrValue), xt);
  //450 (the existing record was not compatible).
  cerr << "increment error: " << rdb->error().name() << "incr value: " << retVal << endl;
  record = rdb->get((char *)&key, sizeOfKey, &sp, NULL);

  // Denormalize a 64-bit number in the network byte order into the native order.
  // (big-endian to little endian)
  currValue = *((uint64_t*)record);

  printf("value after increment (should be 146): %d and size %d\n", currValue, sp);

  rdb->clear();

  /* Their CPP version is broke on little-endian (i.e. Intel) machines
  // the CPP version doesn't work either
  rdb->add("1", "123");
  cerr << "add record error: " << rdb->error().name() << endl;
  rdb->increment("1", 23);
  //450 (the existing record was not compatible).
  cerr << "increment error: " << rdb->error().name() << endl;
  rdb->get("1", &value);
  cout << "key 1 incremented to (should be 146): " << value << endl;

  rdb->clear();
  */

  // test bulk set operation
  map<string,string> recs;
  int64_t key1 = 5;
  int64_t key2 = 6;
  int64_t key3 = 7;
    
  int64_t val1 = 5;
  int64_t val2 = 10;
  int64_t val3 = 15;

  recs.insert( pair<string,string>(string((const char*)&key1, sizeof(int64_t)), string((const char *)&val1, sizeof(int64_t))) );
  recs.insert( pair<string,string>(string((const char*)&key2, sizeof(int64_t)), string((const char *)&val2, sizeof(int64_t))) );
  recs.insert( pair<string,string>(string((const char*)&key3, sizeof(int64_t)), string((const char *)&val3, sizeof(int64_t))) );
  retVal = rdb->set_bulk(recs, xt, true);
  cerr << "retval bulk set: " << retVal << endl;
  cerr << " bulk retval: " << rdb->error().name() << endl;

     // retrieve a record
  record = rdb->get((const char *)&key1, sizeof(int64_t), &sp, NULL);
  cout << "key 5 set to (should be 5): " << *((uint64_t*)record) << endl;

  record = rdb->get((const char *)&key2, sizeof(int64_t), &sp, NULL);
  cout << "key 6 set to (should be 10): " << *((uint64_t*)record) << endl;

  record = rdb->get((const char *)&key3, sizeof(int64_t), &sp, NULL);
  cout << "key 7 set to (should be 15): " << *((uint64_t*)record) << endl;

  value = "";

  vector<string> keys;
  keys.push_back(string((const char*)&key1, sizeof(int64_t)));
  keys.push_back(string((const char*)&key2, sizeof(int64_t)));
  keys.push_back(string((const char*)&key3, sizeof(int64_t)));
  retVal = rdb->remove_bulk(keys,true);
  cerr << " bulk records removed: " << retVal << endl;
  cerr << " bulk retval: " << rdb->error().name() << endl;
  record = rdb->get((const char *)&key3, sizeof(int64_t), &sp, NULL);
  if (record == NULL) {
    cout << "record 7 deleted: success!" << endl;
  } else {
    cout << "record 7 not deleted: failed!" << endl;
  }

  rdb->clear();
  if (!rdb->close()) {
    cerr << "close error: " << rdb->error().name() << endl;
  }

  return 0;
}
